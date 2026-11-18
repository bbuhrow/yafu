# GNFS table rows C150-C160 (Nov 2026)

New `ggnfs_table_Gimarel` rows for 150-160 digits in steps of 2 (`factor/nfs/nfs.c`). They extend the trend of Gimarel's rows (C127-C148) and were chosen by test sieving. The old 150/155/160 rows are removed. The 165+ rows are unchanged and are still old estimates.

| digits | rlim | alim | lpbr/lpba | mfbr/mfba | lambda | siever | start-q | min_rels |
|---|---|---|---|---|---|---|---|---|
| 150 | 12M | 24M | 30/30 | 60/60 | 2.533 | 13e | 6.0M | 78M |
| 152 | 14M | 26M | 30/30 | 60/60 | 2.533 | 13e | 6.5M | 82M |
| 154 | 16M | 28M | 30/30 | 60/60 | 2.533 | 13e | 7.0M | 86M |
| 156 | 18M | 32M | 30/31 | 60/62 | 2.533 | 13e | 8.0M | 114M |
| 158 | 20M | 36M | 30/31 | 60/62 | 2.533 | 14e | 9.0M | 119M |
| 160 | 24M | 40M | 30/31 | 60/62 | 2.533 | 14e | 10.0M | 124M |

q-range (block per siever call) is 5000 for all rows. yafu interpolates lims, start-q and min_rels between rows and takes lpb, mfb, lambda and siever from the nearer row (ties go to the larger row). For example, C155 gets 17M/30M with 30/31 on 13e.

## Method

- **Harness.** `testsieve_params.py` does the same thing as bayestune's `bt_nfs_ref` (factor/tune.c), but for any parameter set:
  - it runs the same `gnfs-lasieve4I1Xe -f q -c R -o out -n i -a job` command as a real job;
  - six threads sieve at once, each one block of about 120 s at its own special-q, spread from q0 to 1.2 × need / yield(q0);
  - siever start-up time is subtracted;
  - sieve time is integrated linearly in q until `need` raw relations.
- **Calibration.** On the local C140 poly with the C140 row, it predicts within about 10% of the real C140 run (0.016 vs 0.0144 thread-s/q; 4.0-5.3 vs 4.9 rels/q average).
- **Machine and polys.** i5-11600K, 6 threads, WSL, avx512 lasieve5 binaries. Polys are the median-Murphy-E NFS@Home polys in `polys/` (C150 E=5.16e-12, C155 2.81e-12, C160 1.51e-12).
- **Times.** All times below are sieve time only, in wall hours at 6 threads. Bayestune multiplies sieve time by 1.43 for a whole job.
- **Noise.** Each block is noisy (yield varies 20-30% between nearby q). Treat differences under about 5% as ties.

### Relations needed (assumed)

The ranking between lpb choices depends on these numbers, and they are the least certain input:

- **lpb 29/29:** 25.2M + 1.0M per digit above C126 (Gimarel trend). That gives 49M at C150, which agrees with the C157 NFS@Home tight bracket (about 52M raw).
- **lpb 30/30:** 74M + 2M per digit above C148 (Gimarel trend; the C140 run filtered on the first try at 58M with 16k excess cycles). That gives 78M at C150 and 98M at C160. The NFS@Home C164 30/30 bracket (75-81M unique, about 90-96M raw) suggests the trend is a little high by C160.
- **lpb 30/31:** 1.27 × the 30/30 figure. This is a guess: Gimarel's own lpb steps are about 1.3× per bit, and NFS@Home's 30/31 jobs at C160-166 used about 105-113M unique.
- **lpb 31/31:** 1.6 × the 30/30 figure. Also a guess.

## Results

### C150 (current yafu row, from `bayes_tune_i5-11600K_6t.txt`: 14e 25M/25M 29/29, 46.4 h)

| siever | lims | lpb | mfb | lambda | need | sieve h | q span |
|---|---|---|---|---|---|---|---|
| 13e | 7.84M/10.8M | 30/30 | 60/60 | 2.533 | 78M | 42.7 | 2.7-69.8M |
| 14e | 7.84M/10.8M | 30/30 | 60/60 | 2.533 | 78M | 53.1 | 2.7-30.7M |
| 13e | 7.84M/16M | 30/30 | 60/60 | 2.533 | 78M | 36.4 | 4.0-56.4M |
| 14e | 7.84M/16M | 30/30 | 60/60 | 2.533 | 78M | 43.4 | 4.0-22.5M |
| 13e | 12M/16M | 30/30 | 60/60 | 2.533 | 78M | 35.1 | 4.0-51.6M |
| 13e | 12M/16M | 30/30 | 60/60 | 2.7 | 78M | 37.2 | 4.0-52.3M |
| 13e | 12M/16M | 30/30 | 60/60 | 2.4 | 78M | 45.9 | 4.0-71.8M |
| 13e | 12M/16M | 30/31 | 60/62 | 2.533 | 99M | 33.7 | 4.0-50.9M |
| **13e** | **12M/24M** | **30/30** | **60/60** | **2.533** | **78M** | **33.6** | **6.0-47.0M** |
| 13e | 16M/24M | 30/30 | 60/60 | 2.533 | 78M | 34.1 | 6.0-47.5M |
| 13e | 24M/24M | 30/30 | 60/60 | 2.533 | 78M | 34.4 | 6.0-44.1M |
| 14e | 16M/24M | 30/30 | 60/60 | 2.533 | 78M | 39.7 | 6.0-23.8M |
| 13e | 7.84M/10.8M | 29/29 | 58/58 | 2.5 | 49.2M | 47.7 | 2.7-77.0M |
| 14e | 7.84M/10.8M | 29/29 | 58/58 | 2.5 | 49.2M | 57.3 | 2.7-32.6M |
| 13e | 7.84M/10.8M | 31/31 | 62/62 | 2.533 | 125M | 49.2 | 2.7-78.7M |
| 14e | 7.84M/16M | 31/31 | 62/62 | 2.533 | 125M | 46.4 | 4.0-24.1M |
| 14e | 16M/24M | 31/31 | 62/62 | 2.533 | 125M | 39.6 | 6.0-23.2M |

### C155

| siever | lims | lpb | mfb | need | sieve h | q span |
|---|---|---|---|---|---|---|
| 13e | 16M/24M | 30/30 | 60/60 | 88M | 70.4 | 6.0-84.9M |
| 13e | 16M/32M | 30/30 | 60/60 | 88M | 64.9 | 8.0-81.9M |
| 13e | 24M/32M | 30/30 | 60/60 | 88M | 63.8 | 8.0-71.4M |
| **13e** | **16M/24M** | **30/31** | **60/62** | **112M** | **60.4** | **6.0-78.3M** |
| 14e | 16M/32M | 30/30 | 60/60 | 88M | 73.8 | 8.0-39.1M |
| 14e | 16M/32M | 30/31 | 60/62 | 112M | 69.2 | 8.0-36.1M |
| 14e | 16M/32M | 31/31 | 62/62 | 141M | 64.0 | 8.0-32.7M |

### C160 (current yafu row, from bayestune: 14e 40M/40M 30/30, need 96.5M, 113.9 h)

| siever | lims | lpb | mfb | need | sieve h | q span |
|---|---|---|---|---|---|---|
| 13e | 16M/32M | 30/30 | 60/60 | 98M | 125.0 | 8.0-149.0M |
| 13e | 24M/40M | 30/30 | 60/60 | 98M | 115.7 | 10.0-122.5M |
| 13e | 24M/32M | 30/31 | 60/62 | 124M | 107.6 | 8.0-126.9M |
| 13e | 24M/40M | 31/31 | 62/62 | 157M | 104.1 | 10.0-117.4M |
| 14e | 24M/40M | 30/30 | 60/60 | 98M | 119.4 | 10.0-61.6M |
| **14e** | **24M/40M** | **30/31** | **60/62** | **124M** | **105.0** | **10.0-53.7M** |
| 14e | 24M/40M | 31/31 | 62/62 | 157M | 108.5 | 10.0-52.6M |

All lambdas at C155/C160 were 2.533.

## Reading

- **Lambda.** 2.533 (Gimarel's value) was the best of the three tried. 2.4 loses about 30%, and 2.7 is slightly worse.
- **Lims.** Gimarel's alim stopped growing at 10.8M from C136 on, which is too small at C150. Doubling alim saves about 20%. Above that the gain is flat: 12M/16M through 24M/24M are tied at C150. The rows grow lims slowly, roughly with the q range a job needs.
- **13e vs 14e.**
  - C150 at lpb 30: 13e is about 15% faster.
  - C155: 13e is still ahead.
  - C160: tied. I chose 14e from C158 because 13e has to sieve q up to 120-150M (3-4× alim). That raises the duplicate rate, which a raw `min_rels` does not account for. 14e covers the job in q 10-54M.
- **lpb.**
  - 30/31 ties or beats 30/30 from C150 up, and is ahead by 5-8% at C155-C160.
  - 31/31 is close at C160 but rests on a weaker guess for relations needed.
  - The rows switch to 30/31 at C156. If full jobs show 30/31 needs more than about 1.35× the 30/30 count, move the switch to C158 or later.
- **Expected gain over the old rows on these polys.**
  - C150: about 28% less sieving (33.6 h vs 46.4 h).
  - C160: about 8% (105 h vs 114 h).
  - The old C160 row (40M/40M, 30/30, 14e) was already close to the optimum. The old C150 row (lpb 29, large lims, 14e) was not.

## Predicted times for full jobs on these polys (6 threads, i5-11600K)

| size | sieve (test-sieve) | whole job (×1.43) |
|---|---|---|
| C150 | 33.6 h | about 48 h |
| C155 | about 60 h | about 86 h |
| C160 | 105 h | about 150 h |

To use the same polys, put only the poly lines (`n`, `skew`, `c0`-`c5`, `Y0`, `Y1`) from `polys/cNNN_*.poly` in `nfs.job`. yafu then fills in the parameters from the table and skips poly selection.

## What the full jobs should report

- Whether filtering succeeds at `min_rels` (78M / 114M / 124M). If filtering fails, yafu sieves more and retries, which gives a bracket. For a tighter bracket, start lower with `-minrels`, for example 0.85×.
- Raw vs unique relations (the duplicate rate on 13e's long q span vs 14e).
- The final q reached, and sieve wall time compared with the predictions above.
