# Bayesian ECM-vs-sieve stopping rule

Status (2026-11-15, branch `ecm_prob`, uncommitted): the core model and a standalone test driver are written. It is not wired into `factor()` yet. That waits on local timing data: NFS runs (`docs/nfs_job_data/local_collection.md`) and ECM curve timings from an extended tune.

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

## 3. Open items

1. **Timing inputs.** Extended (bayes) tune, separate from `tune()`: ECM seconds per curve at each B1 for input sizes at the AVX-ECM limb-count boundaries, at the run's thread count; NFS T(D) from `docs/nfs_job_data` (min_rels x 0.85, test-sieve rate, LA fit, poly select budget); SIQS from the existing tune; per-method thread efficiency.
2. **Integration.** Add `-plan bayesian` (`PRETEST_BAYES`) to `schedule_work`. At each ECM state build options for the standard B1 values with `ecm_plan_curve_params` plus measured costs, and call `ecmbayes_decide` with `ecm_obj.twork`. Run the step it picks (B1 and curves) in place of the next t-level state. Fall back to `normal` without tune data. Log the decision (best B1, cost, pfind, ratio, T(D)). `-work` already feeds `twork`. Add `-bayes_prior random|hard`.
3. **Step choice.** The greedy ranking can pick small B1 values when the quantum is small. That doesn't change the stop much but costs about 1% in E[total] (quantum 10 s against T/1000 at C150). Consider ranking by benefit/cost over a few look-ahead steps, or restricting to B1 >= the current t-level's B1.
4. **Cofactor value.** V(c) ignores further ECM on the cofactor (an upper bound). That is fine in the GNFS range; check near the SIQS crossover.
5. **Prior for real inputs.** Validate on aliquot and Cunningham cofactor corpora: compare the predicted find rate per t-level with the observed one.
