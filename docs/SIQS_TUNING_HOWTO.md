# SIQS parameter tuning: how-to

How to produce new SIQS parameters with `siqs_tune_v2.py`, check them, and
get them into yafu, either per machine (`siqstune()` → `siqs_params.txt`) or
as built-in defaults (`factor/qs/siqs_params.c`).

## The pieces

| piece | what it does |
|---|---|
| `siqs_tune_v2.py` | the optimizer. **Pass 1:** a deterministic grid per digit size. **Pass 2:** a local refine and a check on held-back numbers. |
| `siqs_tune.py` | the older random-search tuner. v2 imports its run cache, output parsing and export, so its `report`/`export` also work on a v2 output directory. |
| `tfsm_compare.py` | one-off experiment: adaptive TFSm vs `-noopt` vs fixed TFSm offsets. |
| `siqstune()` | yafu command. Benchmarks the exported parameter sets on the current machine and writes `siqs_params.txt`. |
| `siqs_params.txt` | per-machine table. When present in the working directory, SIQS takes its parameters from it for every option not given on the command line. |
| `factor/qs/siqs_params.c` | the built-in tables, per pipeline (bucket / direct) and kernel ISA, used when there's no `siqs_params.txt`. |

## Before starting

- **Linux or WSL, run from the repo directory.** Each yafu run gets a pseudo-terminal, which needs Linux.
- **Build yafu first, and don't rebuild during a tuning run.** The run cache doesn't record which build made a run. v2 writes `<outdir>/binary.json` and refuses a different binary, so **use a new `--outdir` for each build.**
- **Keep the machine quiet.** Another job on the CPU skews timings badly: one 4-thread job slowed a whole c70 round by ~25%. For overnight runs, turn off Windows sleep and pause updates. WSL stops when Windows sleeps.
- **Pick one thread count and keep it for the whole of pass 1**, or sizes aren't comparable. Each size's thread count, pipeline and MFBD are stored in `grid_cNN.json`; a later run with different values is refused.
- **Put the output directory outside the repo** (or in `.gitignore`). `runs.jsonl` gets large.
- **Scripts need Unix line endings.** A `.sh` written on Windows passes a stray `\r` argument and every line fails. Fix with `sed -i 's/\r$//' script.sh`.

The examples below use:

```bash
OUT=/mnt/c/Users/bbuhr/projects/factoring/v2_<build-date>
```

## Pass 1: the grid

### What's swept

The grid is 4 × 3 × 3 × 3 × 3 = **324 points per test number**. Values are defined relative to what yafu chooses for that number, so they compare across numbers and sizes:

| param | levels | relative to |
|---|---|---|
| B | ×0.5, 0.65, 0.8, 1.0 | yafu's default B |
| M | 100, 180, 300 | absolute |
| TF | −8, −4, 0 | yafu's default TF |
| TFSm | ≤ c70: F−6, −4, −2 · ≤ c80: F−4, −2, 0 · above: F−2, 0, +2 | F, the formula value (what yafu uses with `-noopt`) |
| NB | ≤ c70: 1, 2, 3 · ≤ c80: 2, 3, 4 · above: 3, 4, 5 | absolute |

Fixed for the whole pass:
- **MFBD = 1.72** (`--mfbd`).
- **The pipeline:** `Mid=sparse,Large=direct,Kernel=cmem` by default (`--pipeline`). It's set explicitly so that changing NB can't switch pipelines, which yafu's automatic choice would do. Below 70 digits MFBD and Large/Kernel have no effect.

Any level list can be overridden with `--B`, `--M`, `--TF`, `--TFSm` or `--NB` (comma lists), e.g. to narrow the grid for larger sizes.

### Step 1: the first number, every point

```bash
python3 siqs_tune_v2.py grid --yafu ./yafu --digits 70 --inputs 1 --threads 1 --outdir $OUT
```

- **Reference first.** Each number starts with a reference run: yafu's defaults with `-noopt` on the fixed pipeline, run twice. It gives the default B, TF and NB, the formula TFSm (F), and the reference time. Every result is a **ratio** to that reference time; below 1.0 means faster.
- **Order and timeouts.** Points run in shuffled order. A run slower than `--timeout-mult` × the reference (default 1.5) is stopped and recorded as "slower than 1.5×".
- **Progress.** Each run prints its ratio, the best so far and an ETA.
- **Resuming.** Interrupting is safe: rerun the same command and it continues from the cache.

### Step 2: look at it, and prune

```bash
python3 siqs_tune_v2.py analyze --digits 70 --outdir $OUT
```

It prints:
- **Per level:** the best, median and worst ratio for each parameter level. A level whose best is far behind the overall best is safe to drop.
- **Pairwise tables:** the best ratio for each cell of B×M, TF×TFSm, B×NB and B×TF. These show interactions, for instance whether the best M changes with B.
- **The best points**, with how many numbers each was measured on.
- **A suggested prune list**, written to `$OUT/prune_c70.txt`. It drops whole levels whose best point is more than 8% slower than the overall best (`--level-tol`), and single points more than 15% slower (`--point-tol`), timed out or failed.

**Edit `prune_c70.txt` before using it.** It's plain text: one `level P=V` or `point B=..,M=..,TF=..,TFSm=..,NB=..` per line, and `#` starts a comment. Delete any line you want to keep running. Noise is ~1–2% per run, so keep the thresholds wide; the cost of keeping a point is one run per further number.

### Step 3: the other two numbers

```bash
python3 siqs_tune_v2.py grid --yafu ./yafu --digits 70 --inputs 2,3 --prune --threads 1 --outdir $OUT
```

`--prune` alone reads `$OUT/prune_c70.txt`. A different file also works, and `{d}` in the path is replaced by the size: `--prune /path/prune_c{d}.txt`. **Use the same level options as step 1**, or the grids won't line up.

Then check the result across all three numbers:

```bash
python3 siqs_tune_v2.py analyze --digits 70 --outdir $OUT --no-write
```

### Step 4: fit

```bash
python3 siqs_tune_v2.py fit --digits 70 --outdir $OUT
```

The fit is a quadratic surface in log B, log M, TF, TFSm and NB, with a separate intercept per number. It uses each number's runs within 25% of that number's best (`--fit-within`), because a quadratic only describes the region near the optimum.

- **Which terms it uses** depends on the levels present: a parameter with one level is left out, two levels give a straight-line term, three or more add a curvature term. Cross terms are kept between every pair of parameters that vary.
- **R² and residual sd:** a residual sd near the run noise (1–2%) means the surface fits well. Much larger means the region is too wide or the surface is lumpy; try `--fit-within 0.15`.
- **Fitted optimum**, relative and absolute for a typical number of that size, with the predicted ratio. The optimum is only searched inside the measured range, never extrapolated.
- **Best measured point** on all numbers, for comparison. If it's far from the fitted optimum, look at the tables before trusting either.
- **Sensitivity:** the predicted ratio one step either side of the optimum for each parameter. A flat direction (≈1.00 both ways) is a parameter that barely matters at this size, a candidate to fix or drop at larger sizes.

The result is saved as `$OUT/fit_c70.json`.

### All the small and medium sizes

A driver for one night, the first number at every size:

```bash
for d in 50 55 60 65 70 75; do
  python3 siqs_tune_v2.py grid --yafu ./yafu --digits $d --inputs 1 --threads 1 --outdir $OUT
done
```

Next day: `analyze` each size, edit the prune lists, then run inputs 2 and 3 (same loop, `--inputs 2,3 --prune`). Then `fit` each size.

Rough cost on an i5-11600K, 1 thread:

| size | run time | input 1 (324 runs) |
|---|---|---|
| c50–c60 | under a second | minutes |
| c65 | a few seconds | ~15 min |
| c70 | ~7 s | ~45 min |
| c75 | ~18 s | ~2 h |

Inputs 2 and 3 cost about the same again, less whatever was pruned. Bad grid corners stop at the 1.5× timeout. c80 on 1 thread is ~4–5 h per number. Use more threads for c80 and up, but then use the same count for every size you'll compare.

## Trends, and larger sizes

```bash
python3 siqs_tune_v2.py trends --outdir $OUT
```

It prints, for each size:
- the fitted optimum (B as a multiple of the default, M, TF and TFSm offsets, NB), its predicted ratio and R²;
- the optimum in absolute terms, and the best measured point.

Then a linear trend per parameter against digits, with predictions for c85/90/95 (`--predict`). Use those to centre a **narrower grid** for the larger sizes, for example with 2–3 levels per parameter around the predicted values:

```bash
python3 siqs_tune_v2.py grid --yafu ./yafu --digits 85 --inputs 1 --threads 6 --B 0.55,0.65,0.75 --M 180,300 --TF -6,-3 --TFSm -2,0,2 --NB 3,4,5 --outdir $OUT
```

Parameters that came out flat in the sensitivity output can be fixed at one level.

## Pass 2: refine and check

```bash
python3 siqs_tune_v2.py refine --yafu ./yafu --digits 70 --threads 1 --outdir $OUT
```

- **Centre.** It starts at the fit's optimum (`--center fit`). Alternatives: the best measured grid point (`--center best`), or explicit absolute values (`--center B=9808,M=190,NB=2,TF=72,TFSm=14`).
- **Search.** It takes `--trials` (40) small random steps, measured on the grid's numbers against their reference times. Sets that come within twice the noise of the best are re-run up to `--max-reps` times. The box allows B ×0.8–1.25, M ×0.7–1.4, NB ±1, TF ±3 and TFSm ±3. At ≥ c70, MFBD is also free within ±0.1. `--free-pipeline` frees Mid/Large/Kernel too.
- **Check on held-back numbers.** Finally it runs `--holdout` (3) new numbers × `--final-reps` (2) rounds, comparing:
  - yafu with no flags (the reference, ratio 1.0);
  - yafu's built-in parameters on the bucket pipeline they were tuned for;
  - the best refined set, plus the most different sets within `--diverse-tol` (5%) of it, `--final-top` (5) in total.

It prints and saves `$OUT/summary_c70.json` (siqs_tune.py's format).

## Export, then siqstune on each machine

```bash
python3 siqs_tune_v2.py export --outdir $OUT
```

This writes `$OUT/siqs_candidates.txt` (and `.json`): for each size and each exported held-back number (`--numbers`, default 1), one row per set. The rows carry the number they were measured on, and that number's partials, fulls, needed, sieve and post-processing times. Rows:

| rank | what | how siqstune runs it |
|---|---|---|
| 0 | yafu with no flags | without flags; reference only |
| 1 | built-in parameters on the bucket pipeline, formula TFSm | forced exactly as listed |
| 2 | the best refined set | forced exactly as listed |
| 3+ | the most different near-best sets (dropped if more than `--max-gap` (10%) slower than rank 2) | forced exactly as listed |

Several output directories can be merged: `--outdir dirA dirB`.

On the machine to tune, put `siqs_candidates.txt` in yafu's working directory and run:

```bash
./yafu "siqstune()" -siqsT 5 -threads 6
```

For each size, every set runs on the same number siqs_tune measured, for `-siqsT` seconds. siqstune extrapolates the partials rate to the partials that set needed, and adds the post-processing time scaled by sieving speed. Small inputs that finish inside the sample are just timed. In tests the estimates for forced sets were 0.4–4.6% off the full runs (one c75 and one c70 table).

- **Rank 0 doesn't compete when extrapolated.** yafu's adaptive TFSm tuning makes its 5-second estimate unreliable (17% low in one test), so it's shown for reference. It only competes when it was timed outright.
- **Output:** the fastest set per size goes to `siqs_params.txt`, which SIQS then uses for any option not given on the command line. Delete the file to go back to the built-in tables.
- **Threads:** siqstune uses the given `-threads`, so tune with the thread count you'll actually run with.

## Updating the built-in tables

siqstune also prints its winners as built-in table rows (bits, B, M, NB, TF, TFSm, MFBD, MFBT; 0 = not tabulated, use the formula or default). They're for the table in `factor/qs/siqs_params.c` that matches the machine's kernel ISA and the rows' pipeline, e.g. `siqs_table_direct_avx512` or `siqs_table_direct`. Paste them in, rebuild, and re-check a few sizes against the old binary.

For a direct-sieve table from the grid alone, without siqstune, `fit_cNN.json` has `optimum_abs`, the absolute values for a typical number of that size.

## Troubleshooting

| symptom | cause / fix |
|---|---|
| `... is not the build this outdir was made with` | the binary changed. Use a new `--outdir`, or `--allow-binary-change` if you really want to mix builds. |
| `--threads ... differs from the existing grid's` | that size was started with different settings. Use them, or a new outdir. |
| `fit: only N runs within 25% of the best, M needed` | too few runs near the best for the terms in use. Run more numbers, or widen `--fit-within`. |
| a point shows `FAILED (log ...)` | yafu didn't factor or crashed. The log is in `$OUT/logs/`. Failed points are excluded and listed in the prune suggestions. |
| `rejected the command line (printed usage)` | the binary is older than one of the options (e.g. `-siqsSPV`, `-siqsMid`). Rebuild. |
| every line of a driver script fails at once | Windows line endings: `sed -i 's/\r$//' script.sh`. |
| `trends: no fit_c*.json` | run `fit` for at least one size first. |

## Files in an output directory

| file | contents |
|---|---|
| `runs.jsonl` | every yafu run: parameters, times, relation counts, effective parameters (cache) |
| `binary.json` | the yafu binary the directory belongs to |
| `grid_cNN.json` | per size: test numbers, reference results, levels run, threads, pipeline, MFBD |
| `prune_cNN.txt` | the editable prune list from `analyze` |
| `fit_cNN.json` | the fitted surface, optimum, sensitivity, best measured point |
| `summary_cNN.json` | `refine`'s result (siqs_tune.py format) |
| `siqs_candidates.txt/.json` | `export`'s output, for siqstune |
| `logs/` | yafu output of failed runs (all runs with `--keep-logs`) |
