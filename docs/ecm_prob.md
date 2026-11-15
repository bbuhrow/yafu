# ECM success probabilities and t-level accounting

Branch `ecm_prob`, November 2026.  Code: `factor/gmp-ecm/ecm_prob.c`,
`include/ecm_prob.h`; accounting in `factor/gmp-ecm/ecm.c`; planner hooks in
`factor/avx-ecm/avxecm_stage2fft.c` (`avxecm_s2fft_plan`, `avxecm_s2fft_info`).

## 1. What changed

Before: `factor()` credited ECM work by looking up "curves per t-level" in
tables in `yafu_ecm.h` that were produced by GMP-ECM for fixed
configurations: GMP-ECM param 0/1/3 with its default B2, and AVX-ECM as
Suyama curves with B2 = 100 B1 and no Brent-Suyama.  None of them described
what AVX-ECM runs by default now (twisted Edwards curves of family (iv), the
FFT stage 2 with a memory-limited automatic B2 and an automatic Dickson
polynomial).  Unknown B1 values were not credited at all, the GPU curves
were recorded with the arguments swapped, and the param 3 table was a copy
of the param 1 table.

Now every batch of curves is recorded with the parameters it actually ran
with:

| run path | family | B2, Brent-Suyama, pairs from |
|---|---|---|
| GMP-ECM (internal or external; sigma given) | Suyama, param 0 | GMP-ECM's stage 2 rules (`ecmprob_gmpecm_stage2`) |
| GPU stage 1 + GMP-ECM stage 2 | param 3 | same |
| AVX-ECM stage 1 + GMP-ECM stage 2 (`-prefer_gmpecm_stg2`) | `-avxecm_curve` | same |
| AVX-ECM | Suyama (`-avxecm_curve 0`) or Edwards (iv) (1, 2) | the run itself (`vec_ecm_main` fills `ecm_obj.last_cp`) |

and adds its expected finds to `ecm_obj.twork[d]` for every factor size
d = 10..90 digits.  From that grid:

* `tlevels[]` (t15 .. t70, as before: 1.00 = one full t-level) are its values
  at multiples of 5;
* the t-level (`total_work`) is the size where the expected finds are 1,
  interpolated in log(work) between whole digits;
* the chance that a d-digit factor was missed is exp(-twork[d]); `factor()`
  prints it after the t-levels, e.g.
  `chance of a missed factor: p30 36.7%, p35 91.8%, p40 99.4%`.

Planning works the same way: `ecm_plan_curve_params(fobj, B1, &cp)` predicts
the parameters a run at B1 would use under the current options (for AVX-ECM
with the automatic B2 it runs the same selection as the stage 2 setup, for
the input's size, threads and memory limit), and `get_curves_for_tlevel` /
`get_curves_required` use the model with them.  `-v -v` prints the plan:

    fac: curves at B1=250000: twisted Edwards family (iv), B2=1.177e+08,
         Brent-Suyama -3, 1.239e+07 pairs, 415.3 curves per t30

The old t-level interpolation (`t(i-1) + 5 * tlevels[i]`, linear in the
fraction of the next level) gave about half a digit more than the work
supports (e.g. t40.6 after exactly one t40 at B1 = 3M); the grid gives
t40.0.  Pretesting targets are therefore reached a little later than before.

## 2. The model

The model is GMP-ECM's (Kruppa, PhD thesis 2010, ch. 5), implemented
independently.  The group order of a curve mod p is taken to be as smooth
as a random integer near x = p e^-delta, delta the curve family's extra
smoothness.  With u = ln x / ln B1:

    P = rho_loc(u)                                            stage 1
      + sum_{B1 < q <= B2} rho_loc(u - ln q / ln B1) / q      stage 2
      + sum_{q > B2} rho_loc(...) / q * E_r[1 - exp(-nr r / q)] Brent-Suyama

rho_loc(u) = rho(u) - gamma rho(u-1) / ln x is the "local" Dickman function
for numbers near x.  The sums are exact over primes below 20000 and
integrals dq / (q ln q) above.  Brent-Suyama: nr is the number of stage 2
pairs (baby steps x giant steps; GMP-ECM's dF^2 k) and r the number of
non-trivial roots of f(a) = +-f(b) mod q, averaged over its distribution:
GMP-ECM's residue classes for GMP-ECM runs, the measured r histograms of
the AVX-ECM FFT stage 2 (`bs_model.py`) for AVX-ECM runs.  A factor of d
digits is taken as p = 10^(d - 0.5), like GMP-ECM.

Dickman rho: exact on [0, 2], then from u rho(u) = int_{u-1}^{u} rho(t) dt
by Simpson's rule with step 1/256 (all terms positive, so the relative
accuracy holds for tiny values): rho(3), rho(10), rho(20) agree with the
reference values to 1e-11, 1e-9, 1e-6.

Two deliberate differences from GMP-ECM 7:

* GMP-ECM evaluates the local correction of the stage 2 integral at x, not
  at the cofactor x / q.  We do the same, since t-levels are GMP-ECM's
  convention (the cofactor version needs ~1.7% more curves).
* For B1 < 20000, GMP-ECM's integral mixes log bases (ln B2 / ln 20000 used
  as a multiple of ln B1), so it counts part of (B1, 20000] twice and skips
  part of (20000, B2].  Its B1 = 2000 and 11000 tables need up to 2x too
  many curves (t15 at B1 = 2000: 34 in GMP-ECM, 25 here).  GMP-ECM also
  cuts rho off at u = 10, which shows as jumps in its far-off-diagonal
  entries (e.g. t65 at B1 = 50000).

## 3. Extra smoothness of the curve families

| family | delta | source |
|---|---|---|
| Suyama (param 0, `-avxecm_curve 0`), param 2 | 3.134 | GMP-ECM `ECM_EXTRA_SMOOTHNESS` |
| param 1 (A = 4d-2, d square) | 3.134 + ln 0.41638 = 2.258 | GMP-ECM (Bouvier, Barbulescu) |
| param 3 (A = 4d-2, d 32 bits; GPU) | 3.134 + ln 0.33048 = 2.027 | GMP-ECM |
| twisted Edwards family (iv) (`-avxecm_curve 1, 2`) | 3.134 + 0.279 = 3.413 | measured, below |

`arith/edwards/fam_delta.py bits n seed` computes #E(F_p) (BSGS in
`edw_ref.py`) for a random curve of each family on n random primes and
delta = sum over l < 200 of (E[v_l(#E)] - 1/(l-1)) ln l.  It is the group
order, not the starting point's order, that matters: the index
#E / ord(P0) is made of small primes that stage 1 covers anyway (on
average ln(index) = 2.1 for Suyama, 2.5 for (iv)), so the part that must
be <= B1, B2 is #E's rough part.

| run | family | E[v2] | E[v3] | delta | delta - Suyama (paired) |
|---|---|---|---|---|---|
| 32-bit p, 40000 | Suyama | 3.335 | 1.685 | 3.143 +- 0.020 | |
| | (iv) | 5.338 | 0.676 | 3.395 +- 0.020 | +0.252 +- 0.028 |
| 40-bit p, 40000 | Suyama | 3.331 | 1.684 | 3.120 +- 0.020 | |
| | (iv) | 5.333 | 0.675 | 3.426 +- 0.020 | +0.306 +- 0.028 |

Suyama's measured delta reproduces GMP-ECM's 3.134 (Kruppa's 2^3.323
3^1.687), so the method is sound; (iv) is anchored at 3.134 + 0.279 (the
mean of the paired differences; E[v2] = 16/3 and E[v3] ~ 0.68 alone give
+0.275).  Family (iv) trades the Suyama curves' 3-torsion for a larger
2-part, net e^0.279 = 1.32x "smaller" orders: 4-6% fewer curves per t-level
than Suyama's at the same B1 and B2 (t30 at 250k with GMP-ECM's stage 2:
407 vs 430).  With the automatic B2 the Edwards curves get a smaller B2
(their stage 1 is cheaper, so stage 2 is balanced against less time) and
need about as many curves per t-level as Suyama's, each one cheaper.

End-to-end check against exact simulation (`edw_ref.py yield`: point
orders of random curves mod random primes, success = order B1-smooth but
for one prime <= B2), model averaged over the same prime range:

| primes | B1, B2 | Suyama model / observed | (iv) model / observed |
|---|---|---|---|
| 40 bits, 60000 | 500, 5e4 | 6.17% / 6.27 +- 0.10% | 6.73% / 6.52 +- 0.10% |
| 48 bits, 30000 | 2000, 2e5 | 4.85% / 4.76 +- 0.13% | 5.21% / 5.29 +- 0.13% |

## 4. Validation against GMP-ECM

`ecmprob_gmpecm_stage2` reproduces GMP-ECM 7.0.6's stage 2 (its bestD with
power-of-2 dF, default k, Dickson degree by stage 2 length) for all twelve
standard B1 values: B2 within 1e-5 relative, and identical dF, k and
polynomial.  GMP-ECM's tables do not depend on the input size (checked for
a C151 and a C99).  Expected curves, model / GMP-ECM -v, params 0, 1, 3, for
35..80 digits: within +-0.5% at B1 >= 250000 wherever GMP-ECM's rho cutoff
does not apply, 1-2% at B1 = 50000; the B1 = 2000 and 11000 columns differ
as described in section 2.

The model with Suyama curves, B2 = 100 B1 and no Brent-Suyama reproduces
the old `avx_ecm_data` table (GMP-ECM -power 1 -param 0 output) to 0.5%:
t30 at 250k 667 vs 671, t40 at 3M 4829 vs 4858, t60 at 260M 119388 vs 119976.

## 5. Tables

Expected curves for a factor of d digits at the standard B1 values; the
diagonal is one full t-level.  Produced with the test driver
(`gcc -O2 -DECMPROB_MAIN -Iinclude factor/gmp-ecm/ecm_prob.c -lm`,
`./a.out std iv gmp` etc.); the automatic-B2 tables with a harness around
`avxecm_s2fft_plan` (IFMA cost model, 8 threads, 256 MB/thread).  The
automatic B2 depends on the input size, threads, memory limit and ISA, which
is why yafu plans per input instead of using a fixed table.

### Edwards (iv), AVX-ECM default: FFT stage 2, automatic B2 and polynomial, 500-bit input

    B1        B2          S
    2e3       1.976e5     -
    11e3      1.784e6     -
    5e4       1.188e7     -
    2.5e5     1.177e8     D3
    1e6       1.298e9     D3
    3e6       5.191e9     D6
    11e6      2.076e10    D6
    43e6      7.787e10    D6
    110e6     1.869e11    D6
    260e6     4.205e11    D6
    850e6     1.345e12    D6
    2900e6    4.877e12    D6

    /*      2e+03    1.1e+04  5e+04    2.5e+05  1e+06    3e+06    1.1e+07  4.3e+07  1.1e+08  2.6e+08  8.5e+08  2.9e+09  */
    /*t15*/ {22, 6, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1},
    /*t20*/ {696, 70, 20, 8, 4, 3, 3, 2, 2, 2, 2, 2},
    /*t25*/ {45794, 1430, 210, 48, 19, 11, 7, 5, 4, 4, 3, 3},
    /*t30*/ {5081965, 45678, 3212, 415, 112, 53, 28, 16, 12, 10, 7, 6},
    /*t35*/ {8.5e+08, 2063428, 66719, 4716, 855, 315, 130, 63, 42, 30, 20, 14},
    /*t40*/ {2e+11, 1.2e+08, 1780465, 67114, 8062, 2284, 735, 288, 169, 110, 66, 41},
    /*t45*/ {6.2e+13, 9.5e+09, 5.9e+07, 1154561, 90547, 19553, 4845, 1517, 780, 456, 238, 133},
    /*t50*/ {2.5e+16, 8.9e+11, 2.3e+09, 2.3e+07, 1182505, 192932, 36429, 9021, 4037, 2107, 959, 471},
    /*t55*/ {1.2e+19, 1e+14, 1.1e+11, 5.5e+08, 1.8e+07, 2155182, 307446, 59726, 23158, 10754, 4243, 1830},
    /*t60*/ {7.4e+21, 1.3e+16, 5.8e+12, 1.5e+10, 3e+08, 2.7e+07, 2874938, 435034, 145554, 59900, 20383, 7693},
    /*t65*/ {5.3e+24, 2e+18, 3.5e+14, 4.3e+11, 5.5e+09, 3.7e+08, 2.9e+07, 3452584, 992866, 360875, 105528, 34696},
    /*t70*/ {4.5e+27, 3.6e+20, 2.4e+16, 1.4e+13, 1.1e+11, 5.6e+09, 3.3e+08, 3e+07, 7294243, 2334993, 584507, 166883},

With 256 MB/thread the automatic B2 at large B1 is ~1800 B1, a third of
GMP-ECM's default: more curves per t-level than GMP-ECM's stage 2 (e.g.
t55 at 110M: 23158 vs 17062 below) but each much cheaper; before this
change those curves were credited as GMP-ECM param 0 curves.

### Edwards (iv), AVX-ECM default with `-avxecm_stg2mem 1024`, 500-bit input

From B1 = 3M on, 256 MB/thread limits the automatic B2.  At 1024 MB/thread
it no longer does: 64 GB/thread chooses the same B2 values, which are where
the predicted time per factor is least (~3200 B1, still below GMP-ECM's
default of ~7000 B1 at 110M).  The B1 <= 1M columns are the same as at
256 MB.

    B1        B2          S
    3e6       7.786e9     D6
    11e6      4.259e10    D6
    43e6      1.384e11    D6
    110e6     3.514e11    D6
    260e6     7.667e11    D6
    850e6     2.439e12    D6
    2900e6    8.796e12    D6

    /*      2e+03    1.1e+04  5e+04    2.5e+05  1e+06    3e+06    1.1e+07  4.3e+07  1.1e+08  2.6e+08  8.5e+08  2.9e+09  */
    /*t15*/ {22, 6, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1},
    /*t20*/ {696, 70, 20, 8, 4, 3, 2, 2, 2, 2, 2, 2},
    /*t25*/ {45794, 1430, 210, 48, 19, 11, 7, 5, 4, 4, 3, 3},
    /*t30*/ {5081965, 45678, 3212, 415, 112, 49, 25, 15, 11, 9, 7, 5},
    /*t35*/ {8.5e+08, 2063428, 66719, 4716, 855, 293, 116, 58, 39, 28, 19, 13},
    /*t40*/ {2e+11, 1.2e+08, 1780465, 67114, 8062, 2115, 651, 264, 154, 101, 61, 38},
    /*t45*/ {6.2e+13, 9.5e+09, 5.9e+07, 1154561, 90547, 18031, 4259, 1379, 705, 416, 218, 122},
    /*t50*/ {2.5e+16, 8.9e+11, 2.3e+09, 2.3e+07, 1182505, 177221, 31816, 8157, 3633, 1913, 876, 433},
    /*t55*/ {1.2e+19, 1e+14, 1.1e+11, 5.5e+08, 1.8e+07, 1972945, 267015, 53776, 20744, 9723, 3859, 1675},
    /*t60*/ {7.4e+21, 1.3e+16, 5.8e+12, 1.5e+10, 3e+08, 2.5e+07, 2484511, 390225, 129854, 53951, 18471, 7015},
    /*t65*/ {5.3e+24, 2e+18, 3.5e+14, 4.3e+11, 5.5e+09, 3.4e+08, 2.5e+07, 3086607, 882599, 323952, 95319, 31536},
    /*t70*/ {4.5e+27, 3.6e+20, 2.4e+16, 1.4e+13, 1.1e+11, 5.1e+09, 2.8e+08, 2.6e+07, 6463291, 2089841, 526431, 151258},

Fewer curves per t-level does not mean less time per t-level: the
automatic B2 minimizes the predicted (stage 1 + stage 2 time) / P(find) at
the factor size the B1 targets (design doc §11.3), so GMP-ECM's larger B2
would take fewer AVX-ECM curves but more total time.  Two caveats: that
choice uses the cost model's time, not measured time, and `s2_prob`'s
simpler probability model (section 7).

### Edwards (iv) with GMP-ECM's default stage 2 (`-prefer_gmpecm_stg2`)

    /*t15*/ {23, 6, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1},
    /*t20*/ {772, 69, 20, 8, 4, 3, 2, 2, 2, 2, 1, 1},
    /*t25*/ {51461, 1410, 201, 47, 19, 11, 7, 5, 4, 3, 3, 2},
    /*t30*/ {5765943, 44958, 3053, 407, 112, 52, 25, 14, 10, 8, 6, 4},
    /*t35*/ {9.7e+08, 2028706, 63083, 4615, 860, 308, 117, 53, 33, 23, 15, 10},
    /*t40*/ {2.3e+11, 1.2e+08, 1676329, 65575, 8102, 2236, 654, 237, 131, 80, 46, 28},
    /*t45*/ {7.2e+13, 9.3e+09, 5.5e+07, 1126798, 90967, 19121, 4280, 1230, 591, 322, 162, 87},
    /*t50*/ {2.9e+16, 8.7e+11, 2.2e+09, 2.3e+07, 1187593, 188470, 31969, 7230, 3014, 1461, 639, 299},
    /*t55*/ {1.4e+19, 9.8e+13, 1e+11, 5.3e+08, 1.8e+07, 2103393, 268282, 47401, 17062, 7330, 2771, 1135},
    /*t60*/ {8.6e+21, 1.3e+16, 5.4e+12, 1.4e+10, 3e+08, 2.6e+07, 2496161, 342295, 105969, 40203, 13082, 4676},
    /*t65*/ {6.2e+24, 2e+18, 3.3e+14, 4.2e+11, 5.5e+09, 3.6e+08, 2.5e+07, 2695744, 715190, 238954, 66687, 20710},
    /*t70*/ {5.3e+27, 3.5e+20, 2.2e+16, 1.4e+13, 1.1e+11, 5.4e+09, 2.8e+08, 2.3e+07, 5204080, 1527383, 364246, 98018},

### Edwards (iv), pairing stage 2, B2 = 100 B1 (`-avxecm_stg2 0`)

    /*t15*/ {21, 7, 4, 3, 2, 2, 2, 1, 1, 1, 1, 1},
    /*t20*/ {693, 79, 25, 11, 7, 5, 4, 3, 3, 2, 2, 2},
    /*t25*/ {45580, 1658, 262, 70, 31, 19, 12, 8, 6, 5, 4, 4},
    /*t30*/ {5056201, 53689, 4101, 630, 196, 95, 48, 27, 19, 15, 11, 8},
    /*t35*/ {8.5e+08, 2451736, 86690, 7410, 1580, 603, 240, 110, 71, 49, 32, 22},
    /*t40*/ {2e+11, 1.5e+08, 2347016, 108485, 15529, 4588, 1421, 527, 297, 187, 108, 66},
    /*t45*/ {6.2e+13, 1.1e+10, 7.8e+07, 1911515, 180551, 40845, 9736, 2878, 1420, 803, 405, 221},
    /*t50*/ {2.5e+16, 1.1e+12, 3.1e+09, 4e+07, 2428564, 416619, 75657, 17669, 7586, 3829, 1683, 810},
    /*t55*/ {1.2e+19, 1.2e+14, 1.5e+11, 9.4e+08, 3.7e+07, 4790169, 656932, 120257, 44701, 20064, 7643, 3231},
    /*t60*/ {7.3e+21, 1.6e+16, 7.9e+12, 2.5e+10, 6.4e+08, 6.1e+07, 6298261, 897231, 287555, 114348, 37561, 13900},
    /*t65*/ {5.3e+24, 2.5e+18, 4.9e+14, 7.7e+11, 1.2e+10, 8.6e+08, 6.6e+07, 7273482, 2002286, 702930, 198368, 63989},
    /*t70*/ {4.4e+27, 4.5e+20, 3.4e+16, 2.6e+13, 2.5e+11, 1.3e+10, 7.5e+08, 6.4e+07, 1.5e+07, 4630328, 1118423, 313349},

### Suyama, GMP-ECM param 0 with its default stage 2 (this model)

    /*t15*/ {25, 6, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1},
    /*t20*/ {848, 74, 21, 8, 5, 3, 2, 2, 2, 2, 1, 1},
    /*t25*/ {57375, 1526, 214, 50, 20, 11, 7, 5, 4, 3, 3, 2},
    /*t30*/ {6499636, 49120, 3274, 430, 118, 54, 26, 14, 10, 8, 6, 4},
    /*t35*/ {1.1e+09, 2233232, 68112, 4909, 906, 323, 122, 55, 34, 23, 15, 10},
    /*t40*/ {2.6e+11, 1.3e+08, 1820033, 70104, 8574, 2351, 684, 246, 135, 83, 47, 28},
    /*t45*/ {8.3e+13, 1e+10, 6e+07, 1209759, 96646, 20179, 4487, 1282, 614, 334, 168, 89},
    /*t50*/ {3.3e+16, 9.7e+11, 2.4e+09, 2.5e+07, 1266039, 199535, 33616, 7558, 3140, 1518, 662, 309},
    /*t55*/ {1.7e+19, 1.1e+14, 1.1e+11, 5.8e+08, 1.9e+07, 2233185, 282851, 49673, 17815, 7630, 2874, 1174},
    /*t60*/ {1e+22, 1.5e+16, 5.9e+12, 1.5e+10, 3.2e+08, 2.8e+07, 2637925, 359495, 110882, 41940, 13597, 4844},
    /*t65*/ {7.3e+24, 2.3e+18, 3.6e+14, 4.6e+11, 5.9e+09, 3.8e+08, 2.7e+07, 2836839, 749783, 249733, 69433, 21489},
    /*t70*/ {6.2e+27, 4e+20, 2.5e+16, 1.5e+13, 1.2e+11, 5.8e+09, 3e+08, 2.4e+07, 5465288, 1598962, 379852, 101859},

(The old GMP-ECM tables in `yafu_ecm.h` came from A. Schindel and the
AVX-ECM table from tbusby; they were removed with this change, the model
reproduces them as described in section 4.)

## 6. Test driver

    gcc -O2 -DECMPROB_MAIN -Iinclude factor/gmp-ecm/ecm_prob.c -lm -o ecm_prob_test
    ecm_prob_test rho                     # rho against reference values
    ecm_prob_test gmp p0 3e6              # GMP-ECM's stage 2 and its -v table (35..80)
    ecm_prob_test avx iv 3e6 5.2e9 -6 5.3e8   # AVX-ECM curves: B1 B2 S pairs
    ecm_prob_test std iv gmp|avx100|edwfft    # table at the standard B1 values
    ecm_prob_test prob iv B1 B2 S pairs bsmodel digits

## 7. Open items

* The AVX-ECM automatic B2 / polynomial choice (`s2_prob` in
  `avxecm_stage2fft.c`) has its own simpler copy of the model (delta = ln 23,
  no local correction, the Suyama delta for Edwards curves too).  It could
  call `ecmprob_prob` with the run's family.
* Family (iv)'s delta could be refined with more samples (+-0.02 now, which
  is +-0.3% in curves per t-level).
* GMP-ECM without the NTT (very large inputs) uses non-power-of-2 dF;
  `ecmprob_gmpecm_stage2` assumes the NTT.
