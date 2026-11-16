# Bayesian ECM-vs-sieve stopping rule

Status (2026-11-15, branch `ecm_prob`, uncommitted):
- Written: the core model with a standalone test driver, the ECM part of the extended tune (`bayestune()`, §3), and `-plan bayesian` in `factor()` (§4).
- Not done: the NFS part of `bayestune()`. Until then, sieve times come from the `tune()` fits in yafu.ini.

Code: `factor/ecm_bayes.c`, `include/ecm_bayes.h`. It depends only on `ecm_prob.h`.

## 1. Model

**Prior.** The prime factors of N in digit band d (10^(d-1) <= p < 10^d, the band that `ecm_obj.twork[d]` describes) are independent Poisson counts with mean

    lambda[d] = ln(d / (d-1))      (Mertens; random integer)

This is truncated below `dmin` digits (trial division, rho and factors already known) and above D/2 digits, since the smallest factor is at most sqrt N. The mass between 90 digits and D/2 is kept in one extra slot. ECM never finds those factors, but they count when conditioning on N being composite. The `HARD` prior puts all the mass within `margin` digits of D/2, for inputs believed to be balanced semiprimes.

**Posterior.** Completed ECM work W[d] (expected finds, i.e. `twork`) thins the process. Undetected factors have mean lambda[d] exp(-W[d]). Conditioning on at least one factor remaining (N is composite), a step that adds work w[d] finds a factor with probability

    P = (1 - e^-Ls) / (1 - e^-L),    L = sum lambda e^-W,   Ls = sum lambda e^-W (1 - e^-w)

As W grows, the remaining mass moves toward D/2, so a deeply pretested input is treated more and more as a hard semiprime.

**Value of a find.** If a d-digit factor is found, the sieve time saved is T(D) - V(D-d), where V(c) = P(c-digit cofactor composite) * T(c). It is averaged over the bands the step reaches. T is a callback: min(SIQS, GNFS) from tune data, or the SNFS time when an SNFS polynomial will be used.

**Rule.** Each failure lowers every later success probability, so this is a monotone stopping problem and one-step look-ahead is optimal. Every candidate B1 is evaluated for a step of equal time (default T(D)/1000, at least `min_curves`). The step with the best benefit/cost runs while that ratio is at least 1; otherwise the sieve starts. Sunk ECM time does not enter the decision.

API: `ecmbayes_prior`, `ecmbayes_option_init` (call once per B1, caches one curve's work grid), `ecmbayes_decide`. `ecmbayes_pfind` and `ecmbayes_saving` are exposed for logging.

## 2. Test driver

    gcc -O2 -DECMBAYES_MAIN -Iinclude factor/ecm_bayes.c factor/gmp-ecm/ecm_prob.c -lm
    ./a.out [-k s_per_1M_curve_at_C150] [-tscale x] [-q quantum_s] [-d from to step] [-v]

It simulates the all-fail path on a model machine:
- GNFS: 3960 s at C120 (the measured local run in `nfs_job_data/c120`), doubling every 7 digits.
- SIQS: doubling every 4 digits below C95.
- ECM curves: Edwards (iv) with GMP-ECM stage 2, k (B1/1e6) ((D+10)/160)^2 s.

It reports the expected total time (ECM plus sieve on N or the cofactor) for the Bayesian stop and for fixed 4/13 and 1/3 targets. Steps are chosen by the same ranking in all three, so the only difference is where they stop. **All the timings are placeholders until tune data exist.**

Default model (k = 0.1, i.e. fast AVX-ECM on 6 threads):

| D | T(D) | Bayes stop | ECM/T | E[total] Bayes | 4/13: stop, E[total] | 1/3: stop, E[total] |
|---|---|---|---|---|---|---|
| 100 | 547 s | t38.2 | 0.29 | 90 | t30.8, 126 | t33.4, 107 |
| 120 | 3960 s | t41.5 | 0.21 | 892 | t37.0, 978 | t40.0, 903 |
| 150 | 7.7e4 s | t47.3 | 0.16 | 16883 | t46.2, 16987 | t50.0, 17619 |
| 180 | 1.5e6 s | t54.1 | 0.13 | 3.04e5 | t55.4, 3.06e5 | t60.0, 3.74e5 |
| 200 | 1.1e7 s | t59.0 | 0.12 | 2.08e6 | t61.5, 2.14e6 | t66.7, 2.97e6 |

What this shows:
- With these placeholder timings, the 4/13 plan is close to optimal from C140 to C180: the Bayesian stop is within about 1 digit and 1% in expected time. The rule gains most at small sizes, where the sieve is cheap relative to ECM setup and 4/13 under-pretests (10-30% in this model), and above C190 (3%).
- Sensitivity matches the planning estimate. Halving T moves the stop about 2 digits (C150: t47.3 to t45.4). ECM 4x slower (k = 0.4) moves it 4-5 digits shallower, and then the gain over 4/13 is 5% at C150 and 17% at C180. So getting the ECM curve cost right matters more than fine NFS calibration.
- The expected total under random-integer inputs is dominated by early finds: P(find) is 0.79 for the first t35-sized step at C150.

## 3. Extended tune: `bayestune()`

This is a separate command (in `factor/tune.c`, listed in `docfile.txt`), so `tune()` is unchanged.

The ECM part runs `ecm_loop()`, the same engine, curve family, stage 2 and `-threads` that `factor()` uses:
- inputs are semiprimes with two factors of half the size, from C60 to C310;
- B1 = 50k, 250k and 1M everywhere, plus 3M for C90–C200;
- each point is 16 curves per thread (two AVX-ECM batches; one batch leaves per-call setup visible at small sizes);
- if a curve finds a factor, the point is redrawn with a new input.

It writes `bayes_tune.txt`:
- a header with the CPU, thread count, date and ECM options;
- one `ecm <digits> <bits> <B1> <curves> <s/curve> <B2>` line per point;
- an `ecm_eff` line: 1-thread against all-thread efficiency at C150, B1 = 250k.

`bayes_tune_load()` and `bayes_ecm_sec_per_curve(fobj, digits, B1)` read the file back:
- between B1 values they interpolate in log-log;
- beyond the measured B1 range they extrapolate with the end segment's slope, clamped to [0.8, 1.5];
- between sizes they interpolate log-linearly, and outside the size range scale as size^2;
- if the thread count differs from the tune, they scale linearly and print a warning.

Measured on an i5-11600K under WSL with 6 threads, 16 curves per thread: 529 s. The file is saved as `docs/nfs_job_data/bayes_tune_i5-11600K_6t.txt`.

Wall seconds per curve with 6 threads:

| C | 50k | 250k | 1M | 3M |
|---|---|---|---|---|
| 60 | 0.0018 | 0.0058 | 0.018 | |
| 100 | 0.0028 | 0.012 | 0.039 | 0.13 |
| 150 | 0.0040 | 0.020 | 0.082 | 0.24 |
| 200 | 0.0057 | 0.033 | 0.138 | 0.41 |
| 310 | 0.015 | 0.059 | 0.257 | |

Compared with one batch per thread (8 curves):
- The per-call setup is visible at small sizes and small B1: one batch read 15-45% slower at C60-C120.
- From C130 up at B1 >= 250k the two runs agree within about ±8%, which is the run-to-run noise. The earlier non-monotone points (C150 vs C160, C190 vs C200) are gone.
- The C150, 1M figure (0.082) is close to the test driver's placeholder (k = 0.1).

**Thread efficiency is 0.60 at 6 threads** (C150, B1 = 250k: 0.064 s/curve on 1 thread, 0.018 on 6), and it was the same with one batch per thread. So it is not setup overhead. Likely causes are AVX-512 all-core clocks and memory bandwidth in the FFT stage 2. This is the size of thread-scaling effect discussed in planning, and it is below the 0.75 that `get_*_time_estimate` assumes for sieving. Measuring all three methods at the run's thread count keeps the ECM/sieve ratio honest.

## 4. `-plan bayesian`

`-plan bayesian` (alias `bayes`) is `PRETEST_BAYES`. It needs `tune()` data (`tune_info` for this CPU) and `bayes_tune.txt`; without either, it prints one notice and follows the normal plan. With both, in `schedule_work()`:
- At each ECM state, `bayes_schedule()` builds one option per standard B1 level that still has curves left:
  - curve parameters from `ecm_plan_curve_params`;
  - seconds per curve from `bayes_ecm_sec_per_curve`;
  - at least 8 x threads curves per step.
- It calls `ecmbayes_decide()` with `ecm_obj.twork`, the random-integer prior, and T(d) = min(SIQS, GNFS) from the tune fits, divided by threads x 0.75 as in `get_*_time_estimate`. The input's own job is timed at the SNFS equivalent size when an SNFS polynomial was chosen.
- It then either runs the chosen level with the chosen curves (capped at that level's maximum) or switches to the sieve.
- The t-level target is disabled while the plan is active, including during the P-1/P+1 states. The normal target is still used only for the trivial-input check.
- The options' work grids are cached per input size. Recomputing them for every level on every step cost about 0.5 s per step.
- **Start-up pause (fixed; it affected every plan).** `init_factor_work()` plans the curves of all 11 ECM levels, and for AVX-ECM each plan runs the FFT stage 2 automatic-B2 search (`avxecm_s2fft_plan` → `s2_auto_B2`). That took 77 ms per level, about 0.84 s on a C91 before trial division started. Three changes in `avxecm_stage2fft.c`, none of which changes any result (704-point grid of size, B1, threads, memory and stage-1 type, identical output):
  - `s2_prob` is evaluated once per (polynomial, t, blocks) instead of once per (k, lean) configuration as well, with the same tie-breaking;
  - `s2_nbaby` is cached per D;
  - `avxecm_s2fft_plan` results are memoized (64 entries).
  A plan now takes 8.5 ms, and repeat calls are free, so the pause is about 0.09 s.
  - The per-level maximum curves are now planned only when a level is first reached: `ecm_level_max()` in autofactor.c, resolved inside `get_next_state()`, and also used by the `-work` preload and `bayes_schedule()`. The pause before trial division is gone, which matters for long lists of small inputs.
  - One side effect: a level's maximum now uses the size of the composite being worked on when the level is reached (for example the cofactor left after rho), not the original input's size.
- Each decision is printed at `-v` and logged to factor.log ("bayes: t.., sieve .. s, best .. curves at B1=.. P(factor) .. ratio ..: run|stop").

The build files changed with it:
- `ecm_prob.c` and `ecm_bayes.c` are now common sources (Makefile `COMMON_SRCS`, VS `factor_common` project). `siqs_demo` links `factor_common`, which calls `ecmprob_init`, so it failed to link before.
- The Makefile link targets now depend on the archive rules instead of re-archiving inline. Under `make -j`, several targets rewriting the same `.a` raced and broke the `ecm_demo` link.
- All four targets build with `make -j12` in WSL. The MSVC and MinGW builds are untested.

First runs (i5-11600K, WSL, 6 threads):

| input | normal plan | bayesian |
|---|---|---|
| C80 = p40 p41 | t24.9, 8.2 s | t28.2, 11.5 s |
| C100 = p30 p70 | 5.4 s (p30 found by ECM) | 8.5 s (p30 found at B1 = 1M) |

How to read these:
- The predicted step costs matched the run: about 0.5 s per 48 curves at B1 = 250k on the C100.
- The C80 stop is deeper than the normal plan because ECM is cheap next to the estimated 7 s SIQS. That is the right choice in expectation over random inputs, but it costs about 3 s on a balanced semiprime.
- Single inputs can't show the expected gain. The comparison needs a set of random composites (or aliquot cofactors) per size, run with both plans.
- At C100 the decision used a SIQS estimate of 432 s, extrapolated from the tune fit. It is worth checking against a real C100 SIQS run.

## 5. Open items

1. **Remaining timing inputs** for `bayestune()`:
   - NFS T(D) from `docs/nfs_job_data`: min_rels x 0.85, test-sieve rate, LA fit, poly select budget;
   - SIQS from the existing tune;
   - per-method thread efficiency.
2. **Integration follow-ups.**
   - Add `-bayes_prior random|hard` (the code always uses the random prior).
   - Benchmark normal against bayesian on sets of random composites per size.
   - Build and test under MSVC and MinGW.
3. **Step choice.** The greedy ranking can pick small B1 values when the quantum is small. That doesn't change the stop much but costs about 1% in E[total] (quantum 10 s against T/1000 at C150). Consider ranking by benefit/cost over a few look-ahead steps, or restricting to B1 >= the current t-level's B1.
4. **Cofactor value.** V(c) ignores further ECM on the cofactor (an upper bound). That is fine in the GNFS range; check near the SIQS crossover.
5. **Prior for real inputs.** Validate on aliquot and Cunningham cofactor corpora: compare the predicted find rate per t-level with the observed one.
