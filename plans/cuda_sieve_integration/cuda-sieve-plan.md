# Plan: cuda-sieve as an external GNFS lattice siever in YAFU

Written 2026-09-21. Progress and decisions live in `status.md`; this file is the
stable reference. Upstream asks are in `cuda-sieve-upstream-requests.md`.

## Goal and scope

Let YAFU's NFS sieving phase call the external CUDA lattice siever
(`kyleaskine/cuda-sieve`, executable `bench`) the same way it calls
`gnfs-lasieve4I1Xe`: one process per range, output concatenated into YAFU's
relation file, progress tracked in `<outputfile>.ranges`.

**Assumptions**
- The cuda-sieve binary is used **unmodified**. Upstream requests are optional
  improvements, not dependencies.
- **No new state-machine state.** CADO-NFS has its own `NFS_STATE_CADO`; this
  work instead adds a backend switch inside the existing lasieve hooks.
- `test_sieve` (poly-select hook) is **out of scope** and assumed disabled or
  reworked separately.
- The non-threadpool sieving path is dead (`USE_THREADPOOL` is hard-defined in
  `nfs_sieving.c`); only the threadpool path is touched.
- Linux first; Windows in the last phase.

## How the two sievers map

| | lasieve (today) | cuda-sieve |
|---|---|---|
| Executable | `<ggnfs_dir>gnfs-lasieve4I<N>e` | user-supplied path to `bench` |
| Range | `-f start -c qrange` = `[start, start+qrange)` | `--qrange S:E`, **inclusive** -> `E = start+qrange-1` |
| Side | `-a` / `-r` | `--sq-side 1` (algebraic) / `0` (rational) |
| Siever version | `I<N>e`, N in 11..16 | `--logI N` (N in 2..20); `-siever N` maps 1:1 |
| Job file | `-a/-r <jobfile>` | `--poly <jobfile>` (reads rlim/alim/lpb/mfb; ignores lambda) |
| Output | `-o file` | `--relations file` (writes `file.part`, renames on completion) |
| Thread id | `-n tid` | `--device D` |
| Cofactoring | inline; optional `-d` batch 3LP + `.raw` | inline on GPU (rho / ECM); no batch step |
| Abort info | `<job>.<host>.last_spq<tid>` | `file.part.ckpt` (`next_q`, `rel_bytes`, ...) |
| Relation format | GGNFS `a,b:rat:alg` | same, special-q included |

Proposed command line:

```
bench --pipeline --cofactor --poly JOB --logI N --qrange S:E --sq-side X \
      --relations rels<T>_<K>.dat --restart --device D  > rels<T>_<K>.dat.log 2>&1
```

## Hook inventory (line numbers approximate, from the provided snapshot)

| File | Function / place | Change |
|---|---|---|
| `cmdOptions.h/.c` | `NUMOPTIONS 131`; `OptionArray`, `OptionHelp`, `needsArg`, `LongOptionAliases`; handler chain indexed `OptionArray[i]`; defaults near line 1442 | Append options 131, 132; bump to 133; extend all four arrays in lock-step; handlers; defaults. Names <= 19 chars, args <= 255 |
| `driver.c` | option -> `nfs_obj` copy (~1804-1811, ~1956-1959) | Copy the two new options |
| `factor.h` | `nfs_obj_t` | Add `cuda_sieve[GSTR_MAXSIZE]`, `cuda_dev[GSTR_MAXSIZE]` |
| `nfs_impl.h` | `nfs_job_t`, `nfs_threaddata_t` | `siever_kind`; per-thread device index, log name, completed-q delta |
| `nfs.c` | `check_for_sievers` (called at ~301) | cuda branch; fix uninitialized `found` and overlapping `sprintf(name, "%s.exe", name)` |
| `nfs.c` | `get_ggnfs_params` sets `sievername` (~2393); `copy_job` (~2536) | Route through one helper |
| `nfs_filemanip.c` | `parse_job_file` sets `sievername` from `# siever:` (~1440) | Same helper |
| `nfs_sieving.c` | `lasieve_launcher` (~2691) | Backend branch: command, log redirect, exit decode, salvage |
| `nfs_sieving.c` | `nfs_sieve_start` (~394) | Skip 3LP batch init and `gpu_device_init` for cuda; range sizing; device map |
| `nfs_sieving.c` | `nfs_sieve_sync` (~667) | Abort/partial handling from `.ckpt`; fix `rels%d.dat` name |
| `nfs_sieving.c` | `nfs_afb_prime_cache` (~2018) | Early return for cuda (P1); optional roots1 cache (P4) |

Not touched: `nfs_sieve_dispatch`, `.ranges` format, `savefile_concat`, msieve
filtering, CADO glue.

## Options (proposed)

- `-cuda_sieve <path>`: path to the `bench` executable; presence selects the cuda backend.
- `-cuda_dev <list>`: comma-separated device indices, e.g. `0,1`. Default `0`.
- Existing `-siever N` selects `--logI N`; `-use_gpudev` stays ECM-only.
- Selecting cuda forces `nfs_batch_3lp` off, with a printed notice.

## Phases

### P0: plumbing, no behavior change without the option
- Add options and fields as in the inventory.
- One helper `nfs_set_sievername()` used by `get_ggnfs_params`, `parse_job_file`
  and `copy_job`.
- `check_for_sievers`: with cuda selected, check the binary exists and skip the
  lasieve scan (otherwise yafu falls back to SIQS on a cuda-only machine).
- Fix the pre-existing `check_for_sievers` bugs and the `rels%d.dat` name.

**Done when:** yafu builds; lasieve command lines are byte-identical before and
after (`-v3` prints `syscmd`); `-h` lists the new options.

### P1: sieving on one GPU
- `lasieve_launcher`: cuda branch builds the command with `snprintf` and fails
  loudly on truncation (buffer is 1024). Skip `-d` and `process_batch`.
- Redirect siever output to `<outfile>.log`; print its tail on nonzero exit.
- Decode `system()` status portably (`WEXITSTATUS` on POSIX).
- `nfs_sieve_start`: for cuda set `is_3lp = 0`, skip relation-batch init and
  `gpu_device_init`.
- `nfs_afb_prime_cache`: return early for cuda.
- Pre-check the job file: duplicate recognized keys and missing `skew` are
  fatal in cuda-sieve's parser.
- THREADS=1, device from `-cuda_dev`, GNFS with algebraic special-q.

**Done when:** a c100-c120 factors end to end; `bench --check-relations` passes
on a sample output; msieve filtering accepts the relations; relations/q is in
line with a lasieve run on the same job.

### P2: stops, failures, range bookkeeping
- Completion test: final output file exists.
- **Stop detection:** a clean stop exits 0 but leaves only `.part`. `system()`
  hides Ctrl-C from yafu, so detect "exit 0, no output file, `.part` present"
  and set `NFS_ABORT`.
- **Salvage:** parse `.part.ckpt` (`next_q`, `rel_bytes`, `relations`), truncate
  `.part` to `rel_bytes`, rename to the output name, and record
  `next_q - start` (capped at qrange) as the delta `nfs_sieve_sync` writes to
  `.ranges`. Replaces the `.last_spq` search for this backend. No checkpoint
  means discard.
- Exit codes: 3 fatal with a "rebuild with wider BN_LIMBS" message; 4 failed
  range (not logged, retried by later dispatch); 5 keep salvaged relations and
  warn that mfb/bucket sizing is off; anything else fatal.
- Delete leftover `.part`, `.ckpt`, `.lock`, `.log` for each output name before
  launch. Always pass `--restart`.
- SNFS: `--sq-side 0` from `poly->side`.

**Done when:** Ctrl-C mid-range yields a correct partial file and `.ranges`
line, and resuming yafu continues from the right q; `kill -9` of the siever
loses at most one flush; a forced exit 3/4/5 is handled as above.

### P3: range sizing and multi-GPU
- Size `thread_qrange` from measured seconds per q (`test_time` already exists)
  so a launch runs about 5-10 minutes (tunable). This replaces the 3LP
  `MIN(1000, ...)` and 1M-raw-rels logic for this backend.
- Thread `i` uses `--device list[i % n]`; warn when THREADS exceeds the device
  count (VRAM, host contention).
- Note the `timeout` check only fires between ranges, so it gets coarser.
- Decide `--blocking-sync` by measurement.

**Done when:** two-GPU run scales; range durations sit near target; timeout
behavior documented.

### P4 (optional): factor-base cache
- Measure in-process GPU factor-base generation per launch first. Proceed only
  if it is a meaningful fraction of a launch.
- Run `fbgen_gpu` once per job, gated by `keep_afb`, pass `--fb1`. YAFU must
  validate the roots1 header (`# Roots for polynomial`, `# lim`, `# maxbits`)
  itself, since cuda-sieve appears to check only `maxbits`.

### P5: validation and packaging
- Yield/duplicate comparison against lasieve on the same jobs (GNFS and SNFS).
- Full factorization on several sizes; msieve filtering and sqrt.
- Windows build and `system()` quoting; `--stop-file` note (TerminateProcess
  cannot checkpoint).
- Docs, help text, pin a tested cuda-sieve commit.

## Risks

- **Host CPU contention** costs cuda-sieve up to ~20% throughput (its STATUS
  findings 53/96). Avoid CPU lasieve alongside by default.
- **Scale/allowance derive from each launch's first q**; many short launches
  differ slightly from one long band. Measure the yield effect.
- **Duplicates** from the full special-q-side factor base (15-25% per its STATUS)
  are expected; msieve filtering removes them.
- **Multiple processes per GPU** compete for VRAM (several GB each) and host.
- **User job files** with duplicate keys or no `skew` fail in cuda-sieve.
- **Stale roots1 cache** would silently give wrong roots (P4 only).
- **Requirements:** NVIDIA sm_80 or newer, `lpb <= 64`, `mfb <= 128`, at most
  three large primes per side, norms within the build's `BN_LIMBS`.

## Out of scope
- `test_sieve`, poly-select coupling, `lasieve_launcher_tdata` (unreferenced in
  the provided files).
- Any change to cuda-sieve itself.
- Replacing YAFU's own GPU batch cofactoring for lasieve.
