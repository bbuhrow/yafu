# RNS modular multiplication / squaring: design notes

Status: decisions recorded (§0); scalar reference in arith/rns/.

## 0. Decisions (2026-11-06)

* Q1 layout: in the small regime, 8 independent numbers per vector (e.g. 8
  curves, same N), with channels as the loop index, so each extension-matrix
  load is reused 8 times. In the Yao regime, 8 moduli per vector (one number
  per call). Caveat: the Yao channel count is then fixed at a multiple of 8,
  so we can't always sit at the optimum n of the paper's Fig. 2.
* Q2 second extension: Shenoy-Kumaresan redundant channel, with
  m_r = 2^52 (arithmetic is just the low 52 bits; IFMA madd52lo). The first
  extension stays inexact (Bajard): Q^ = Q + alpha*M_B with alpha <= n_B, absorbed by the bounds.
* Q3 size: initial use is 512-2048-bit ECM. The long-term aim is FFT/NTT
  channel products in the Yao regime: still ~2n^2 extension work, but
  asymptotically fast per-channel products, and smaller transforms keep more
  bits of precision.
* Q4 ISA: AVX-512 (IFMA, plus F-only FP emulation) and AVX2.
* Q5 API: a standalone RNS API (arith/rns/rns.h).

Scalar reference (arith/rns/, `make` in WSL, `./rns_test [iters] [seed]`):
* Moduli are FCFS 2^52 - d (d odd), pairwise coprime but not prime. N must
  be coprime to M_B*M_C. For ECM, a shared factor is a found factor.
  d_max is 557 at 2048 bits and 1187 at 4096 bits, so 2^52 = d folds stay cheap.
* Channel counts (nB/nC): 512b 11/10, 1024b 20/20, 2048b 40/40, k = nB+2.
* Residues are partially reduced (< 2^52), as IFMA needs. One lazy lo/hi
  accumulator pair per output channel, reduced once (red116).
* Verified against GMP: random mul/sqr/add/sub chains, worst-case inputs
  h*k*N - 1, long squaring chains, and an output bound < k*N check, at 64-4096 bits.

Gandino reorganization applied (gandino_precomputation_rns.pdf, §2.2):
* -N^{-1} folded into the first extension's B_i^{-1} scaling (f1_b).
* N * M_B^{-1} * C_j^{-1} folded into the first extension matrix (ctx->BtoCN,
  per N, O(nB*nC) setup). The extension sum and the T * C_j M_B^{-1} term
  then share one lazy accumulator per C channel. That saves a product and a
  reduction per channel, and Q^ is never formed.
* C-base inputs are pre-scaled by C_j^{-1}, so the second extension needs no
  diagonal scaling.

AVX-512 IFMA kernel (arith/rns/rns_avx512.c), 8 numbers per vector:
* Bit-identical to the scalar reference in all lanes. Row blocking RB = 4
  (8 accumulators), table entries via broadcast.
* Rocket Lake i5-11600K, gcc 11 -O3, sqrmod ns per number:

  | bits | scalar ref | AVX-512 IFMA | GMP mul+tdiv_r |
  |------|-----------:|-------------:|---------------:|
  |  512 |        281 |         28.2 |            109 |
  | 1024 |        732 |         73.4 |            330 |
  | 2048 |       2773 |        234   |           1021 |
  | 4096 |       9854 |        879   |           3169 |

  At 2048 bits this is about 6600 IFMA per 8 numbers, running at roughly
  0.75 of the 1/cycle 512-bit IFMA peak. -O2 is ~1.6x slower (no unroll of
  the RB loop). RB = 2..8 are within noise of each other.
* The next comparison that matters is AVX-ECM vecsqrmod52 at the same sizes
  (also 8 lanes, CIOS ~2n^2 IFMA pairs). The op counts are nearly equal, so the
  RNS case rests on its lack of carry chains and its N-independent structure.

Benchmark vs AVX-ECM and GMP (arith/rns/rns_bench.c, `make bench`), same
machine. Times are per number (vector times / 8). powm uses a full-size
exponent shared by all lanes and the same sliding-window rule for both vector
engines (AVX-ECM's vecslidingmodexp). Every lane is checked against mpz_powm.

| bits | RNS ch | AE core | sqr RNS | sqr AE | mul RNS | mul AE | powm RNS | powm AE | mpz_powm | GMP/RNS | AE/RNS |
|-----:|-------:|---------|--------:|-------:|--------:|-------:|---------:|--------:|---------:|--------:|-------:|
|  256 |   6+6  | fixed416  |  14.0 |  14.3 |  15.4 |  14.1 |    4.6 us |    3.4 us |    7.5 us | 1.65 | 0.75 |
|  512 |  11+10 | fixed624  |  28.1 |  19.8 |  29.0 |  25.0 |   18.0 us |   12.5 us |   45.4 us | 2.52 | 0.69 |
| 1024 |  20+20 | fixed1040 |  74.1 |  48.3 |  73.0 |  46.0 |   94.3 us |   57.1 us |  271.7 us | 2.88 | 0.61 |
| 1536 |  30+30 | bfips     | 146.2 | 110.7 | 144.8 | 145.6 |  262.8 us |  215.5 us |  862.0 us | 3.28 | 0.82 |
| 2048 |  40+40 | bfips     | 234.9 | 166.1 | 237.0 | 200.1 |  576.1 us |  408.4 us | 1926.1 us | 3.34 | 0.71 |
| 3072 |  60+60 | bfips     | 500.3 | 340.9 | 484.6 | 432.1 | 1735.6 us | 1252.5 us | 5783.5 us | 3.33 | 0.72 |
| 4096 |  80+79 | bfips     | 816.4 | 582.8 | 852.3 | 745.8 | 3926.8 us | 2847.4 us |12659.7 us | 3.22 | 0.73 |

(sqr/mul in ns.) Reading:
* Both vector engines beat mpz_powm: RNS by ~3.3x at 1.5-4k bits, AVX-ECM by
  ~4.5x. AVX-ECM is ~1.4x faster than the RNS small regime.
* The IFMA counts are equal: CIOS is ~4n^2 per 8 numbers, and the two RNS
  extensions are 2*(n+1)*n lo/hi pairs, also ~4n^2. The gap is (a) the O(n)
  RNS overhead (3n channel mulreds plus 2n red116 reductions, roughly 30% on
  top of the matrix work at n = 40) and (b) squaring: AVX-ECM's sqr is ~15-20%
  cheaper than its mul, while RNS gets nothing from squaring in this regime.
  Against mul-vs-mul the gap is only ~0.85-0.9x at >= 1.5k bits.
* The small regime is therefore at best a tie with CIOS on the CPU. The
  expected win is in the Yao regime. Rough count at 2048 bits, n = 8 channels
  per base of 5 words, all in one vector: 5 multi-word channel modmuls (~300
  IFMA) plus the two extensions with single-word B~ constants (~160 IFMA) plus
  reductions, ~600 IFMA per number versus ~800 for CIOS. Squaring also helps
  there, since the channel products are multi-word.

## 0b. Yao regime design (8 + 8 multi-word channels)

Layout: one number per vector call. Base B is 8 channels b_i = 2^w - db_i,
one per SIMD lane, each held as s words of 52 bits (w = 52s); base C is the
same; plus the redundant channel r = 2^52, a scalar. Storage is word-major,
v[(base*s + word)*8 + lane], followed by r. Residues are kept lazily
reduced: any value < 2^w congruent mod the channel modulus.

Capacity: 8w >= maxbits + 2*log2(h*k) + 1 with k = nB + 2 = 10 and h = 8,
so maxbits <= 416s - 14: s = 2: 818, 3: 1234, 4: 1650, 5: 2066, 6: 2482,
8: 3314, 10: 4146 bits. The granularity is coarse (416 bits per word of
channel width). w < 52s is possible later at the cost of shifted folds.

Moduli: plain FCFS 2^w - d, d odd ascending. B gets the 8 smallest d, C the
next 8 (a split, not Yao's interleave), so every constant has a uniform sign:
* first extension: |B_i|_{c_j} ~ B~_ij = prod_{k!=i} (dc_j - db_k) > 0
* second extension: |C_j|_{b_i} ~ -C~_ji, C~_ji = prod_{k!=j} (dc_k - db_i) > 0
Measured (scratch script, w = 104..520): dmax ~ 63, v <= 40.5 bits,
v' <= 39.1 bits (interleaving gives ~36/34 but mixed signs). Every extension
constant is one 52-bit word, so each extension product is one IFMA lo/hi
pair per word, with no masked or signed accumulation.

mulmod (per number; ⊙ = lane-wise multi-word product, fold = reduction mod 2^w - d):
1. B: xi = fold(fold(x_B ⊙ y_B) ⊙ f1),  f1_i = |-N^{-1} B_i^{-1}|_{b_i}
2. first extension: Q_j = fold(sum_i xi_i * B~_ij). This is a cross-lane
   broadcast of xi_i (vpermq per word) times a constant column vector,
   accumulated lazily into s+1 words.
3. C (pre-scaled): xi'_j = fold(fold(x_C ⊙ y_C) * g1_j + Q_j * g2_j), with
   g1 = |C_j M_B^{-1}|, g2 = |N M_B^{-1} C_j^{-1}|. Both products share one
   lazy accumulator. (Gandino's matrix folding is not used here: it would
   make the B~ constants full-width, which loses the Yao advantage.)
4. r: U_r = T_r M_B^{-1} + sum_i lo52(xi_i) |B_i|_r N M_B^{-1} (mod 2^52).
   alpha' = (sum_j lo52(xi'_j) |C_j|_r - U_r) M_C^{-1} (mod 2^52), where lo52
   of a radix-2^52 residue is just its low word.
5. second extension: S_i = sum_j xi'_j * C~_ji, then
   U_i = fold(K_i - S_i + alpha' * |-M_C|_{b_i}), with K_i a precomputed
   multiple of b_i >= max S (about 2^(w + v' + 3)).
Per mulmod: 5 full channel products (2 in B, 3 in C), as in Yao's 5n count,
plus 2 x 8 x s IFMA lo/hi pairs for the extensions. Squaring helps in the
x ⊙ x steps (2 of the 5).

Scalar reference (rns_yao.h, rns_yao_ref.c, rns_yao_test.c; GMP per
channel): passes mul/sqr/add/sub chains, worst-case inputs h*k*N - 1,
squaring chains and the output bound < k*N, from 200 bits up to the s = 10
capacity limit (4146 bits). With plain FCFS, db = 1..~33 and dc = ~33..~87,
v = 40-44 and v' = 39-42 bits (FCFS+ in the script gave ~1-3 bits less).
Setup checks every Bt/Ct constant against the true CRT constant mod the
channel modulus.

If more channels push the constants past 52 bits, Yao §3 "Further
Bit-length Reduction" applies. All moduli are odd, so every difference
(dc - db) is even, and a product of n-1 of them has at least n-1 trailing
zero bits (usually more). Store B~ >> z (the truncated form) and restore the
shift on the product: accumulate per shift class and shift once, or fold the
shift into the constant's word alignment. Not needed at 8 + 8 channels
(v <= 44).

AVX-512 IFMA kernel (rns_yao_avx512.c): one number per call. Channel
products are fully unrolled column sums (lo/hi IFMA into 64-bit columns,
normalized once), specialized per s = 1..10 through an always-inline
template with constant S. Each fold is three branch-free rounds (provably
enough; extra rounds are no-ops). The extensions use a memory broadcast of
one channel's digit times a column of Yao constants. AVX-ECM's vecmul52_n
was not reused: it pads to 4-word blocks (64 products instead of 25 at
s = 5) and returns normalized output, while the C step needs two products
summed before normalization. Bit-identical to the reference at every size;
powm matches mpz_powm.

Benchmark (rns_bench, same machine, per number; powm with a full-size
exponent, the same sliding window for every engine):

| bits | s | sqr small | sqr Yao | sqr AE | mul small | mul Yao | mul AE | powm small | powm Yao | powm AE | mpz_powm | Yao vs AE | Yao vs GMP |
|-----:|--:|----:|----:|----:|----:|----:|----:|-------:|-------:|-------:|--------:|----:|----:|
|  256 | 1 |  14.5 |  60.2 |  12.2 |  14.2 |  59.4 |  12.5 |   4.6 us |  18.6 us |   3.4 us |     7.4 us | 0.18 | 0.40 |
|  512 | 2 |  28.6 |  88.0 |  19.1 |  28.2 |  86.6 |  21.8 |  18.2 us |  53.7 us |  12.5 us |    45.5 us | 0.23 | 0.85 |
|  800 | 2 |  54.7 |  86.1 |  32.2 |  51.4 |  83.3 |  36.2 |  50.7 us |  82.3 us |  31.7 us |   145.9 us | 0.39 | 1.77 |
| 1024 | 3 |  73.5 | 153.9 |  48.8 |  73.8 | 137.0 |  46.5 |  90.8 us | 144.1 us |  58.8 us |   265.5 us | 0.41 | 1.84 |
| 1536 | 4 | 144.2 | 164.8 | 115.9 | 144.0 | 140.3 | 141.9 | 262.0 us | 249.4 us | 214.1 us |   829.2 us | 0.86 | 3.32 |
| 2048 | 5 | 233.6 | 166.7 | 166.4 | 236.0 | 184.3 | 199.4 | 570.1 us | 411.8 us | 406.8 us |  1902.7 us | 0.99 | 4.62 |
| 3072 | 8 | 484.9 | 308.8 | 340.6 | 482.7 | 334.7 | 427.1 |1738.9 us |1129.5 us |1255.6 us |  5896.2 us | 1.11 | 5.22 |
| 4096 |10 | 830.3 | 440.1 | 581.0 | 819.7 | 483.4 | 741.7 |3864.8 us |1984.3 us |2835.2 us | 12630.6 us | 1.43 | 6.37 |

(sqr/mul in ns.) Reading:
* The crossover Yao > small regime is at ~1.5k bits. Yao vs AVX-ECM breaks
  even at 2048 and wins 1.43x at 4096; vs mpz_powm it is 6.4x at 4096.
* Below ~1.5k bits Yao has a ~60-85 ns floor: per-call overhead of 5 folds
  (3 rounds each), 2 extensions, the redundant channel's horizontal sums and
  the scalar alpha' path, with only 8 lanes of work per call. One number per
  call cannot amortize it.
* The apparent Yao "sqr slower than mul at s = 3, 4" was a measurement
  artifact: the first timed loop for each new s ran on cold code and tables.
  rns_bench now warms up first.

Optimization pass 1 (still bit-identical to the reference):
1. Fold fast path. After fold round 1 the remaining top is tiny, so rounds 2
   and 3 matter only if digit 1 carries past 2^52. One compare-mask plus a
   predicted branch replaces two full carry chains per fold. Instrumented
   (-DRNSY_COUNT_SLOW): the slow path fired 0 times in 20k-call loops at every
   size. This was the biggest win (2048 sqr 167 -> 125 ns).
2. Rotated constant tables (s >= 2). F_i = f * 2^(52i) mod m for
   f1, g1, g2, so a product by a per-N constant is sum_i a_i F_i: s+1
   columns, with no high half and no fold round 1. Cost 3*s^2 vectors per
   ctx (19 KB at s = 10). The reference uses the same formula. Gain 5-9%.
3. Squaring accumulates off-diagonal terms directly in acc, then doubles
   (no separate off[] array).
Tried and reverted: issuing the C products alongside the B chain (no gain,
more register pressure at s = 10). Two independent chains interleaved
across calls gain only ~20% at s = 2 and ~0 at s >= 3. The kernel is
throughput-bound: at s = 5 sqr, ~430 IFMA (port 0 only on client cores)
and ~450 add/and/shift (mostly port 5) in ~520 cycles.

After pass 1 (rns_bench, per number):

| bits | s | sqr Yao | sqr AE | mul Yao | mul AE | powm Yao | powm AE | mpz_powm | Yao vs AE | Yao vs GMP |
|-----:|--:|------:|------:|------:|------:|---------:|---------:|----------:|----:|----:|
|  512 | 2 |  60.9 |  22.5 |  61.9 |  22.3 |   39.7 us |   12.5 us |    45.2 us | 0.31 | 1.14 |
|  800 | 2 |  61.4 |  34.3 |  62.5 |  36.9 |   61.0 us |   31.9 us |   147.1 us | 0.52 | 2.41 |
| 1024 | 3 |  80.7 |  46.2 |  80.9 |  46.7 |  100.0 us |   57.0 us |   274.4 us | 0.57 | 2.74 |
| 1536 | 4 |  99.1 | 113.2 | 104.4 | 133.0 |  176.8 us |  205.3 us |   835.7 us | 1.16 | 4.73 |
| 2048 | 5 | 115.8 | 163.1 | 124.9 | 200.2 |  291.8 us |  407.8 us |  1906.4 us | 1.40 | 6.53 |
| 3072 | 8 | 218.8 | 339.5 | 248.8 | 427.6 |  816.8 us | 1252.1 us |  5827.5 us | 1.53 | 7.13 |
| 4096 |10 | 303.6 | 600.8 | 369.3 | 747.7 | 1499.7 us | 2856.3 us | 12728.6 us | 1.90 | 8.49 |

Column-level product fold (RNSY_COLUMN_FOLD in rns_yao_avx512.c, default
1; -DRNSY_COLUMN_FOLD=0 restores the exact kernel). The two variable
products (x*y in B and C) are folded at the column boundary. Low columns
pass through unnormalized; each high column c folds as lo52(c)*d (IFMA
lo/hi) plus (c >> 52)*d one digit up, so one carry chain over s+1 columns
replaces the 2s-column chain. That trades ~4s add/and/shift for s IFMA per
fold. The residues are still valid (congruent, < 2^w, normalized digits,
same bounds) but differ from the reference. rns_yao_test then checks
"same value mod N, every channel encodes the same integer (B CRT vs
pre-scaled C channels vs r), digits normalized, value < bound*N", with a
self-test that corrupted C/r channels and bound violations are rejected.
rnsy_vec_exact() reports the mode.

Measured gain: small, as expected for a kernel that is nearly IFMA-bound.
Steady-state sqr/mul at s = 2..5 is within ~1-3% (2048: sqr 117 -> 116,
mul 128 -> 124 ns); at s = 10 sqr 308-320 -> 306, mul 351-360 -> 344 ns
(~3-4%). powm vs AVX-ECM: 1.43x at 2048, 1.88x at 4096 (noise-level
change). Most remaining normalization feeds IFMA inputs, which must be
normalized digits, so the larger structural lever is IFMA count, not carry
chains.

Crossover vs AVX-ECM is now between 1024 and 1536 bits. Remaining ideas:
* Below ~1.2k bits the ~60 ns floor is per-call overhead. A multi-number
  kernel (2-4 independent numbers interleaved) would recover ~20% at s = 2.
* Normalization is ~40% of the instructions. Dropping exact bit-agreement
  with the reference would allow column-level folds (no full 2s-column
  carry chain), with the tests switched to congruence plus bound checks.
* Merge the two horizontal sums (U_r and alpha') into one reduction tree.
* Karatsuba for the two variable products at s >= 8.
* w < 52s for finer capacity steps (3072 bits pays for s = 8 now).
* Granularity: 3072 bits needs s = 8 (capacity 3314), so ~8% is wasted;
  w < 52s would help.

Fold: X -> (X mod 2^w) + (X >> w) * d, repeated until X < 2^w. It is
deterministic, so the scalar reference (GMP) and the vector kernel produce
identical residues.

## 0c. Two Yao kernels: one number per call vs 8 numbers per call

Both sit on one setup (rns_yao.h / rns_yao_ref.c), generalized to n channels
per base (rnsy_fixed_init_n(maxbits, n); rnsy_fixed_init(maxbits) = n = 8).
* rnsy_v*_avx512 (rns_yao_avx512.c): one number per call, n = 8,
  channels in lanes. Lowest latency, for when one curve must go fast.
* rnsyb_v*_avx512 (rns_yaob_avx512.c): 8 independent numbers mod the same
  N per call (one per lane), channels as the loop index, any n.
  rnsyb_choose_n(maxbits) picks n. The layout is
  v[((base*n + ch)*s + t)*8 + lane], with r at v[2ns*8 + lane], and
  rnsyb_vset_lane / vget_lane convert to and from the single-number form.
Shared helpers (normalization, folds, unrolled products, RNSY_COLUMN_FOLD)
are in rns_yao_vec.h.

Generalizing n uncovered a sign rule the n = 8 case had hidden:
|C_j|_{b_i} == (-1)^(n-1) prod_{k!=j} (dc_k - db_i). For odd n the second
extension adds the sum directly (fx->cneg = 0) instead of K - sum.

Constant sizes (FCFS, split assignment): single-word without shifting up
to n = 9 (v ~ 51 bits). With Yao's trailing-zero trick, applied as one
common shift per output channel (Bt[i][j] = |B_i|~ >> zB[j], with the
summed extension shifted back by zB[j] before folding), n = 10 fits at
every s and n = 11 at most s. Setup returns NULL when nothing fits.

Batch-kernel specifics:
* Every modulus and constant is a broadcast (an embedded {1to8} operand).
  The redundant-channel sums are vertical accumulations, with no horizontal
  reductions, and alpha' is a vector.
* alpha' * |-M_C|_{b_i} comes from a 16-entry table per digit through
  vpermi2q (requires n + 2 <= 16), replacing 2s IFMA with s shuffles on the
  non-IFMA port.
* Tested bit-exact against the general-n reference (RNSY_COLUMN_FOLD = 0)
  and by consistency checks (default) for n = 2..11 at 200-4096 bits,
  including shifted-constant setups and odd n. Every lane's powm is checked
  against mpz_powm.

n sweep (batch sqr ns/number, Rocket Lake; noise ~ +-5%): the best n is
3 (512 bits), 4 (800), 5 (1024), 5-6 (1536), 8 (2048), 9-10 (3072),
9-10 (4096). The cost model n(5s^2 + 4s + 14) + 4n^2 s was fit to this
sweep and picks within ~5% of the best at every size.

rns_bench (per number; Yao1 = one-number kernel, Yao8 = batch at the cost-model n/s):

| bits | Yao8 n/s | sqr Yao1 | sqr Yao8 | sqr AE | powm Yao1 | powm Yao8 | powm AE | mpz_powm | AE/Yao1 | AE/Yao8 | GMP/Yao8 |
|-----:|----:|------:|------:|------:|----------:|----------:|----------:|-----------:|-----:|-----:|-----:|
|  512 | 3/4 |  62.2 |  24.6 |  22.8 |   42.0 us |   16.2 us |   13.0 us |    47.8 us | 0.31 | 0.80 | 2.95 |
|  800 | 4/4 |  68.2 |  36.7 |  35.3 |   62.6 us |   35.2 us |   31.4 us |   146.5 us | 0.50 | 0.89 | 4.17 |
| 1024 | 5/4 |  77.8 |  43.9 |  46.6 |   99.1 us |   55.2 us |   56.7 us |   259.8 us | 0.57 | 1.03 | 4.71 |
| 1536 | 6/5 |  96.4 |  76.6 | 109.5 |  171.1 us |  142.9 us |  204.7 us |   825.8 us | 1.20 | 1.43 | 5.78 |
| 2048 | 8/5 | 120.6 | 116.2 | 166.6 |  290.0 us |  279.0 us |  403.7 us |  1902.6 us | 1.39 | 1.45 | 6.82 |
| 3072 |10/6 | 233.8 | 197.6 | 339.3 |  802.8 us |  724.0 us | 1243.1 us |  6024.2 us | 1.55 | 1.72 | 8.32 |
| 4096 |10/8 | 304.8 | 292.4 | 587.7 | 1507.5 us | 1400.6 us | 2860.2 us | 12667.8 us | 1.90 | 2.04 | 9.04 |

Reading:
* Batching removes most of the one-number kernel's small-size floor:
  2.5x faster at 512 bits, 1.8x at 1024. The crossover against AVX-ECM
  moves from ~1.3k bits down to ~1k bits.
* At 2048 bits with the same n = 8, the two kernels are equal (both are
  IFMA-bound). The batch gains at 3-4k bits come from choosing n (10
  instead of 8, with smaller s) and from finer capacity steps.
* Below ~1k bits AVX-ECM's fixed-size cores still lead (0.8-0.9x).
* At s = 10 with n = 8, the batch kernel is slower than the one-number
  kernel (~350 vs ~305 ns sqr). Its working set (x, y, z, scratch and
  tables ~ 60 KB) exceeds the 48 KB L1. The cost model avoids that
  configuration.

Notes from the Yao paper (article-2394.pdf):
* Algorithm 3 (from Gandino) stores the C-base residues pre-scaled,
  X'_j = |X * C_j^{-1}|_{c_j}. The second-extension input is then the C
  output itself, so no separate diagonal scaling step is needed. We do the same.
* Yao uses Kawamura alpha estimation for both extensions; we use it for neither
  (the first extension is inexact, the second uses the redundant channel).
* Complexity: 2n^2 (v x w) + 5n (w x w) channel multiplies. v (the bit-length of
  the B~ constants) fits v ~ 0.13n^2 + 4.2n + 7.3 for FCFS+ moduli, so v <= 52
  holds up to about n = 8 (v ~ 49). Target: a fast `x*y mod N` / `x^2 mod N`
engine for ECM, P-1, P+1 on N from roughly 1k bits up to very large sizes.

## 1. Overview (as proposed)

1. Purpose: a modmul/modsqr backend for ECM / P-1 / P+1 stage 1 and stage 2.
2. Algorithm: RNS Montgomery multiplication, Bajard-style two-base
   (B, B') with base extensions, plus Gandino et al.'s constant
   folding / precomputation.
3. Small regime: K channels of 52-bit moduli; all channel arithmetic
   with AVX-512 IFMA, or IFMA emulated with SIMD FP (the existing AVX-ECM
   non-IFMA path).
4. Base extension in the small regime uses lazy reduction: sum
   the K products in 64-bit lanes, letting carries pile up in the spare high
   bits, and reduce once at the end.
5. Large regime (Yao, Fan, Cheung, Verbauwhede, IEEE TC 2014): a few
   multi-word channels m_i = 2^x - c_i (c_i small, pairwise coprime).
   Channel arithmetic reuses the AVX-ECM multi-precision kernels.
   The crossover point is to be determined.
6. Setup is split into a fixed part (depends only on the bases, i.e.
   on M, M' for a given size class) and a variable part (depends on N).

## 2. Cost model (why a crossover exists)

Let n = size of N in 52-bit words.

Small regime, K ~ n channels per base. Each modmul is:

* channel products: 2K mults (B and B'), and
* two base extensions, each a K x K matrix-vector product: about 2K^2 MACs.

Total is about 2n^2 MACs, the same count as schoolbook CIOS. The win is not in
the op count. It is that there are no carry chains, every MAC is independent,
reduction is lazy, and the regularity suits SIMD and GPU well.

Large regime, K channels of w = n/K words each:

* channel products: about 2K * w^2 = 2n^2/K (less with Karatsuba/FFT), and
* base extensions with Yao's small constants: 2 * K^2 (w x 1) products, about
  2Kn word MACs.

Total is about 2n^2/K + 2Kn. That is minimized near K ~ sqrt(n), giving
O(n^1.5). With subquadratic channel multiplication the extension term
dominates and K should shrink. This is the basis for both crossovers: single
precision to Yao, and later schoolbook/Karatsuba channels to FFT/NTT channels.

Memory note for the small regime: each extension matrix is K^2 words (K = 200
gives 320 KB per matrix, K = 1000 gives 8 MB). Past L2 size a single
modmul is bandwidth-bound unless matrix loads are amortized over several
independent multiplies (see Q1). The Yao regime avoids this, because its
matrices are K x K of small scalars with K small.

## 3. Yao constants: where the "N x 1" comes from

With m_j = 2^x - c_j and M_i = prod_{k != i} m_k:

    M_i mod m_j = prod_{k != i} (c_j - c_k)  mod m_j        (for j in the other base)

So the extension constants are products of K-1 small differences. They stay
single-word as long as roughly (K-1) * log2(max |c_j - c_k|) < 52. This
couples K to the spread of the c's. The base B' needs its own
c'_j interleaved with the c_i. The user's 2^x - y form also gives cheap
channel reduction (fold hi*c + lo), and cheap Kawamura alpha estimation
(top bits of xi_i / 2^x).

To confirm from the paper: how Yao handles the per-channel scaling
xi_i = x_i * |M_i^{-1}|_{m_i}. In general that is a full w x w modmul per
channel (K of them), which is the same order as the main channel products.
Gandino-style folding can merge it with the -N^{-1} step, so there is one
full product per channel, not two.

## 4. Setup split (item 6), refined

Fixed tables (per base pair B, B'; reusable for every N of that size class):

* moduli, reduction constants (c_i, Barrett/Montgomery per channel if not
  pseudo-Mersenne)
* |M_i^{-1}|_{m_i}, |M_i|_{m'_j}, |M'_j^{-1}|_{m'_j}, |M'_j|_{m_i}
* |M|_{m'_j}, |M'|_{m_i}, |M^{-1}|_{m'_j} for alpha correction / scaling
* the redundant-channel or alpha-estimation constants

Variable tables (per N), mostly Gandino folding:

* |-N^{-1} * M_i^{-1}|_{m_i}: one constant per channel, so the first
  extension's diagonal scaling and the -N^{-1} multiply become one product.
* |N|_{m'_j} folded with |M^{-1}|_{m'_j}, and the second extension's
  diagonal scaling |M'_j^{-1}|_{m'_j} folded into |N * M^{-1} ...|
  per channel.
* Montgomery-domain constants: |M mod N| and |M^2 mod N| in both bases, for
  conversions in and out.

Open: Gandino also folds diagonal scalings into the K x K matrices.
That makes the matrices N-dependent, costing O(K^2) per N, but it removes
O(K) work per modmul. For ECM, where one N runs millions of modmuls,
O(K^2) per N is negligible. So the split may be "fixed = moduli plus raw
tables; variable = folded matrices". Measure before deciding.

## 5. Lazy reduction bounds (item 4)

A 52-bit digit in a 64-bit lane leaves 12 spare bits, not 13. With IFMA
each product splits into lo52 and hi52, accumulated in separate 64-bit
accumulators:

* lo accumulator: sum of K values < 2^52, which holds up to K < 2^12 = 4096.
* hi accumulator: the same.
* final value = hi*2^52 + lo (~116 bits). Reduce it mod m_j once, cheap if
  m_j = 2^52 - c_j (fold twice).

So the small regime is safe up to K = 4096 channels (about 200k bits), far
past the expected Yao crossover. With FP emulation the same layout holds,
since the existing AVX-ECM FP path already produces lo/hi 52-bit halves.

Suggestion: use pseudo-Mersenne 52-bit channels 2^52 - c_j in the small
regime as well (Yao's idea applied at word level). Channel modmul then needs
no Barrett/Montgomery constant. There are thousands of usable pairwise-coprime
c_j < 2^20.

## 6. Correctness bounds that drive K

* Represented values live in [0, cN) for a small c, with no final
  subtraction. Need M > c^2 N * (something) and M' > c N + K*m-ish,
  depending on the extension variant. Pin down the exact inequalities for the
  chosen variant (Bajard first extension inexact by up to K*M, second exact).
* ECM/P+1 also need modadd/modsub in RNS. These are channel-wise and
  cheap, but each add widens the represented bound (sub adds a multiple of N).
  Either reserve headroom in M for a fixed number of adds between multiplies,
  or provide a "re-tighten" op. This affects how many extra channels are needed.
* Output conversion (CRT back to binary) is only needed for gcd / final
  results, so it can be slow. Mixed-radix or a plain CRT with GMP is fine.

## 7. Open design questions

Q1. Batch layout. Which do we vectorize across?
  (a) Channel-major: one number per call; SIMD lanes hold 8 channels. The
      matrix-vector product vectorizes over output channels, with a broadcast
      xi_i.
  (b) Number-major: 8 independent numbers (8 ECM curves, same N) per
      vector, with channels as the loop index and moduli broadcast. The base
      extension becomes a K x K by K x 8 matrix product (GEMM-like), so every
      matrix load is reused 8 times, which fixes the bandwidth problem in §2.
      This matches AVX-ECM's existing 8-curves-per-vector model.
  For Yao channels, (a) is the natural fit for the existing vec_bignum
  kernels: 8 lanes = 8 different moduli 2^x - c_i, the same shape as AVX-ECM's
  per-lane moduli. Neither existing reduction path does a per-lane c
  (isMersenne is a single c for all lanes).

Q2. Which second-extension method? Shenoy-Kumaresan redundant channel, or
  Kawamura FP alpha estimate? With near-2^52 or near-2^x moduli, alpha
  estimation is almost free and saves one channel.

Q3. Size range and the third regime: hundreds of millions of bits, with
  FFT/NTT per channel. At that point one large FFT multiply plus Montgomery
  (gwnum/IBDWT) is strong, so the RNS case rests on parallelism (GPU, many
  cores) and on not needing a special-form N. Is the
  goal general N at that size (e.g. ECM on large cofactors), or is this mainly
  a GPU path?

Q4. ISA scope: AVX-512 IFMA plus AVX-512F FP emulation to start. Also
  AVX2 (4 lanes, FP-emulated) and NEON, in line with the SIQS multi-ISA work?

Q5. Consumers / API: plug into AVX-ECM as an alternative to
  vecmulmod52 at larger sizes (so it keeps the vec_bignum interface), or
  provide a standalone `rns_ctx` API with its own ECM/P-1 drivers?

Q6. Squaring: in the small regime, RNS squaring costs the same as
  multiplication (each channel is a single product). It only saves in the Yao
  regime, where each channel's product is multi-word. So the squaring advantage
  appears only above the crossover.

## 8. Proposed first steps (pending answers)

1. A scalar reference implementation (uint64 + `__int128`, GMP for checking):
   small-regime Bajard + Kawamura with Gandino folding, plus a random modmul
   test against mpz. This pins down the bounds in §6.
2. A scalar reference for the Yao regime, using GMP mpn for channel products.
3. AVX-512 IFMA small-regime kernel in the chosen layout (Q1), then FP
   emulation.
4. Vector Yao kernel built on AVX-ECM multi-precision routines, plus a
   per-lane-c pseudo-Mersenne reduction.
5. A crossover benchmark against AVX-ECM vecmulmod52 and GMP mpn.
