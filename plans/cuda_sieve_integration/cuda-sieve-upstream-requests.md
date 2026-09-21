# cuda-sieve: requests from the YAFU integration

YAFU will drive `bench` (and optionally `fbgen_gpu`) as an external process, the
same way it drives `gnfs-lasieve4I1Xe`. The first integration targets the
**unmodified** binary; nothing below blocks it. These requests would remove
workarounds and make the integration more robust. Priority: **H** = removes a
correctness/robustness risk, **M** = removes a workaround, **L** = polish.

Status of what was read: README, RUNBOOK, STATUS, `bench_main.cu` (CLI),
`ckpt.h`, `poly.c`, `fb_cado.c`, `fbgen.c` at commit `2802631`. Not run on a GPU
yet; statements about runtime behavior come from the docs and code comments.

## How YAFU will call it

```
bench --pipeline --cofactor --poly JOB.job --logI N \
      --qrange START:END --sq-side {0|1} \
      --relations rels<T>_<K>.dat --restart --device D [--blocking-sync]
```

- YAFU hands out half-open ranges `[start, start+range)`; it passes
  `END = start+range-1` because `--qrange` is inclusive.
- One process per GPU; YAFU maps its worker threads to `--device`.
- Many launches per factorization (minutes each), each with its own output file.
- YAFU concatenates finished output files into its own relation file and logs
  `side,start,delta,relations` per range so a restart can skip finished work.

## Interface assumptions we rely on (please keep stable, or tell us)

1. **Clean stop exits 0 and leaves `NAME.part` + `NAME.part.ckpt`; `NAME` appears
   only when the band completes.** YAFU distinguishes "finished" from "stopped"
   by whether `NAME` exists.
2. **`NAME.part.ckpt` is text `key = value`** with at least `next_q`, `rel_bytes`,
   `relations`; `next_q` is the first special-q not yet accounted for and
   `rel_bytes` is the valid prefix of `.part`. YAFU will truncate `.part` to
   `rel_bytes` and take `next_q - start` as the completed delta.
3. **Exit codes:** 0 ok/stopped, 1 failed, 3 build too narrow (`BN_LIMBS`),
   4 watchdog stall, 5 degraded (relations kept in `.part`, not committed).
4. **Unknown `.job` keys are ignored; duplicate recognized keys are fatal.**
   YAFU job files will be checked against this.
5. **Relation lines are GGNFS/msieve format, special-q included**, so they can
   be appended to YAFU's msieve savefile unchanged.
6. **`.roots1` header** (`# lim`, `# maxbits`, polynomial comment line) is stable.

## Requests

### R1 (H): validate a `--fb1` file against the job
As far as I could find, `--fb1` loading checks `maxbits` (and prime powers
against it) but not that the header's polynomial and `lim` match the current
`--poly`/`alim`. A stale cache would silently yield wrong-side roots. Please
reject a mismatch (or add `--fb1-check`). YAFU will validate the header itself
in the meantime, but a fixed-format check in the binary is safer.

### R2 (H): unambiguous "stopped" signal
Exit 0 for a clean stop is deliberate (BOINC), but a wrapper must infer the stop
from the missing output file. Please add an opt-in, e.g. `--stop-exit N`
(default unchanged), or a documented final line on stdout such as
`RESULT status=stopped next_q=... relations=...`.

### R3 (M): machine-readable end-of-band summary
`--summary FILE` (or a final `RESULT` line) with: `status`
(`complete|stopped|degraded|failed`), `next_q`, `nq_done`, `relations`,
`rel_bytes`, `ms_per_q`, `wall`. Removes `.ckpt` parsing for the normal path and
gives YAFU the s/q it needs to size ranges (see R8) and test-sieve scores.

### R4 (M): `--version`
Print commit, `BN_LIMBS`, `CF_LMAX`, `PIPE_K`, `GPU_ARCH`, CUDA runtime; exit 0
with no GPU access. Lets YAFU pin a minimum version and refuse a build that is
too narrow for the job before launching hours of work.

### R5 (M): `--list-devices`
Print device index, name, compute capability, total/free VRAM; non-zero exit if
none. YAFU needs a default device list and a preflight; today it can only find
out by failing the first launch.

### R6 (M): `--plan` (dry run)
Parse the job, derive geometry, slab plan and estimated VRAM, print them, and
exit before allocating. Would give YAFU a VRAM preflight per job without
running `testsieve.sh`.

### R7 (M): commit partial output on stop
`--commit-on-stop`: after a clean stop, truncate `.part` to `rel_bytes` and
rename it to `NAME`. YAFU would then treat a stopped band like a short finished
one and only need `next_q` (R3) for its range bookkeeping.

### R8 (L): startup cost visibility
Every launch regenerates the algebraic factor base on the GPU unless `--fb1` is
given. Please print (and include in R3) FB-generation seconds separately from
sieving seconds, so YAFU can size ranges to amortize startup. A managed cache
(`--fb-cache-dir`, keyed by polynomial hash, lim, maxbits, with atomic writes)
would remove the need for YAFU to run `fbgen_gpu` and manage multi-GB files.

### R9 (L): quiet mode
`--quiet` to suppress the periodic progress line. With several workers it
interleaves on the console; YAFU will redirect to a per-thread log meanwhile.

### R10 (L): binary name and install target
`bench` is a very generic executable name. A `cuda-sieve` name (or
`make install` producing `cuda-sieve` and `cuda-sieve-fbgen`) would make
discovery in a directory less ambiguous.

### R11 (L): document `--qrange` boundaries
Confirm in RUNBOOK that `MIN:MAX` is inclusive on both ends, that every root of
every prime q in range is processed, and that adjacent ranges
`[a, b]` and `[b+1, c]` neither overlap nor skip. YAFU relies on disjoint
coverage for its `.ranges` bookkeeping.

## Observations (no action requested)

- The band-wide `scale`/`allowance` are derived from the first q of each launch
  (`params from q=...`). YAFU's many short launches will each derive their own;
  yield differences from this should be small but we will measure them.
- Duplicate relations from the full special-q-side factor base (15-25% per
  STATUS) are expected; msieve filtering removes them.
- Host CPU contention costs measurable throughput (STATUS finding 53/96), so
  YAFU will avoid running CPU sievers next to it by default.
