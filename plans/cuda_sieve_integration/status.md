# Status: cuda-sieve integration into YAFU

Bridge document between phases and conversations. **Update it at the end of every
session** (see "How to use"). Plan: `cuda-sieve-plan.md`. Upstream asks:
`cuda-sieve-upstream-requests.md`.

Last updated: 2026-09-21 (P0 + P1 patch built and runs; not yet run against cuda-sieve itself)

## Start here for Phase 2

The P0/P1 patch is applied, and the yafu binary builds and runs with it.
**Not yet done: any run of the actual cuda-sieve `bench` binary.** So Phase 2
work (salvage from `.ckpt`, exit-code policy, `.ranges` bookkeeping) has two
parts in sequence:

1. First confirm P1 end-to-end on a GPU machine with cuda-sieve built: a
   c100-c120 job, `bench --check-relations` on the output, msieve filtering,
   and a yield comparison against a lasieve run on the same job. This exercises
   the P0/P1 code for the first time against the real binary, not the mock
   script the patch was tested with — expect the command line, exit codes, or
   `.ckpt` format to need small corrections against reality.
2. Once that runs clean, do the Ctrl-C-mid-range test (checks stop detection,
   already in the patch) and then build the P2 salvage/retry logic per the
   plan's P2 section and the "Open items" below.

Bring `cuda-sieve-p0-p1.patch` (already applied — for reference) and this file;
pull the patched `nfs_sieving.c`/`nfs.c`/`nfs_impl.h` from wherever the patch
was applied, since this conversation does not have them.

## How to use

- **Start of a session:** attach or paste this file plus the plan; name the phase
  you are working on. Also attach the source files listed under "Files needed"
  for that phase; line numbers in the plan drift.
- **End of a session:** update the phase table, "Decisions", "Verified / not
  verified", "Open items", and add a "Handoff" entry.
- Keep entries factual. Mark anything not run on real hardware as unverified.

## Phase status

| Phase | Status | Notes |
|---|---|---|
| P0 plumbing | **built** | `cuda-sieve-p0-p1.patch` applied; yafu builds and runs with it |
| P1 one-GPU sieving | **built, unverified against cuda-sieve** | code path runs; **never invoked the real `bench` binary, only the mock script** |
| P2 stops / failures / ranges | partly in P1 | stop detection and cleanup are in; **salvage of partial output, exit-code policy (retry/keep) still open**; blocked on a real cuda-sieve run first (see "Start here") |
| P3 range sizing, multi-GPU | not started | |
| P4 factor-base cache (optional) | not started | measure first |
| P5 validation / Windows / docs | not started | |

## Decisions made

1. Use the **unmodified** cuda-sieve binary (`bench`). Upstream requests are
   optional.
2. **Backend switch inside the lasieve hooks**, not a new state like
   `NFS_STATE_CADO`.
3. **`test_sieve` is out of scope** (assumed disabled or reworked with poly
   select).
4. Only the **threadpool** path is touched (non-threadpool path is dead).
5. Options (proposed, not yet implemented): `-cuda_sieve <path>` selects the
   backend; `-cuda_dev 0,1` lists devices; `-siever N` -> `--logI N`.
6. Selecting cuda forces `nfs_batch_3lp` off (cuda-sieve cofactors inline).
7. Ranges map `[start, start+qrange)` -> `--qrange start:start+qrange-1`
   (cuda-sieve is inclusive). Compute the end in 64-bit.
8. Always pass `--restart`; yafu owns range bookkeeping.

**Deviations of the P0/P1 patch from the plan** (2026-09-21)
9. No `siever_kind` field in `nfs_job_t` and no `copy_job` change: the backend
   is the predicate `NFS_USE_CUDA(fobj)` (`nfs_obj.cuda_sieve[0] != '\0'`).
10. `parse_job_file` is untouched: every sieving path calls `get_ggnfs_params`
    afterwards, and `sievername` is now set only by `nfs_set_sievername()`
    (used by `get_ggnfs_params`). The filter-retry call site does not sieve.
11. The job-file pre-check was dropped. A duplicate key or missing `skew` makes
    cuda-sieve exit 1 at once and the tail of `<out>.log` is printed (tested
    with the "duplicate rlim field" message).
12. Pulled from P2 into P1: a clean stop exits 0 with only `<out>.part`; the
    launcher detects that, sets `NFS_ABORT`, and discards the partial output.
    Without it yafu would relaunch the next range after every Ctrl-C.
13. In P1 every nonzero exit aborts (as for lasieve). Exit 3, 4 and 5 only
    print a specific message.
14. `check_for_sievers` with cuda selected checks only that binary and does
    **not** fall back to SIQS if it is missing.
15. Also fixed in `check_for_sievers`: uninitialized `found`, and the
    overlapping `sprintf(name, "%s.exe", name)` (now `strcat`).

## Verified / not verified

**Verified by reading code or docs**
- cuda-sieve commit `2802631` (2026-09-16); executable `bench`; CLI per its
  usage text; output is GGNFS/msieve format with special-q included.
- Clean stop exits 0, leaves `NAME.part` and `NAME.part.ckpt`, no final `NAME`.
  Exit codes: 3 build too narrow, 4 watchdog stall, 5 degraded.
- `.ckpt` is text `key = value` with `next_q`, `rel_bytes`, `relations`,
  `nq_done`, fingerprint, scale/allowance.
- cuda-sieve's job parser ignores unknown keys, requires `skew`, and treats
  duplicate recognized keys as fatal. Tested with the CPU-side parser
  (`fbgen`) on a synthetic job using the keys yafu writes.
- YAFU side: `USE_THREADPOOL` defined at `nfs_sieving.c:44`;
  `lasieve_launcher_tdata` is not referenced in the provided files;
  `check_for_sievers` would revert to SIQS without lasieve binaries;
  `sievername` is set in `get_ggnfs_params` and `parse_job_file`;
  `fill_job_file` appends only missing keys; the sieve side is not stored in the
  job file (`job->poly->side`).
- Pre-existing YAFU bugs seen: `nfs_sieve_sync` abort fallback opens
  `rels%d.dat` (tpool names are `rels%d_%d.dat`); `check_for_sievers` reads an
  uninitialized `found` and uses overlapping `sprintf`.

**Verified for the P0/P1 patch (2026-09-21, Linux, no GPU)**
- `patch -p1` and `git apply --check` succeed on the CRLF originals; the result
  is byte-identical to the working copy.
- Option tables: `OptionArray`, `OptionHelp`, `needsArg`, `LongOptionAliases`
  all have 133 entries; handler indices go to 132; names are <= 18 chars.
- `cuda-sieve-p0-p1-tests.tar.gz` (`sh tests/run_tests.sh <patched dir>`)
  passes: the launcher against a mock `bench` (command line, `[q0, q0+qrange-1]`
  mapping incl. 64-bit end, device round robin, success, clean stop, exit 1/3/5,
  silent exit 0, killed, missing binary, stale files, too-long command line,
  bad logI); `nfs_set_sievername` gives the same lasieve name as before
  (also with `-DWIN32`); `check_for_sievers` cuda and lasieve branches; the
  `-cuda_dev`/`-cuda_sieve` handlers.

**Verified since (2026-09-21)**
- Patch applied to the real tree; yafu builds and the resulting binary runs.
  (Which platform/compiler, and whether `-cuda_sieve`/`-cuda_dev` were
  exercised on the CLI, wasn't reported — confirm at the start of the next
  session if it matters.)

**Not verified**
- The cuda-sieve `bench` binary has **never actually been invoked** by this
  code — P0/P1 testing used a shell-script mock. The real command line, exit
  codes, log format, and `.part`/`.ckpt` naming are unconfirmed against the
  real binary.
- Nothing has been run on a GPU; runtime behavior of cuda-sieve is from docs
  and code comments only.
- That a completed band removes its `.ckpt` (per `pipeline.cuh`, not observed).
- That `--fb1` loading does not validate polynomial/lim (only `maxbits` seen).
- Real yafu job files (only the writer code was read).
- Yield, speed, VRAM and startup cost of cuda-sieve on the target hardware.
- Behavior with multiple processes per GPU.

## Files needed per phase

- **P0:** `cmdOptions.h`, `cmdOptions.c`, `driver.c`, `factor.h`, `nfs_impl.h`,
  `nfs.c`, `nfs_filemanip.c`, `nfs_sieving.c`.
- **P2-P3:** the **patched** `nfs_sieving.c`, `nfs_impl.h`, `nfs.c` (apply
  `cuda-sieve-p0-p1.patch` first, or attach the patched files), plus the
  cuda-sieve RUNBOOK sections "Stopping and resuming" and the CLI usage.
- **P4:** the above plus `bench/FBGEN_GPU.md`, `bench/fbgen_gpu.cu`.

Project files currently available: `nfs_sieving.c`, `nfs.c`, `nfs_impl.h`,
`factor.h`, `nfs_filemanip.c`, `cmdOptions.h`, `cmdOptions.c`, `driver.c`.

## Open items

- Upstream requests not yet sent to the cuda-sieve author (file is ready).
- Option names `cuda_sieve` / `cuda_dev` are used in the patch; confirm or rename.
- Known P1 limitations: `thread_qrange` is still yafu's default (launches may be
  too short for a GPU: startup overhead, P3); `test_sieve` is not handled and
  would break with the cuda backend (`-b -k -c 0 -F`, `strstr` on the name);
  `-v` verbosity is not forwarded; paths with spaces are not quoted (same as
  lasieve); siever output goes to `<out>.log`, deleted on success, kept on
  failure.
- Default mapping `siever N` -> `--logI N`; a GPU may prefer larger logI, tune
  in P3/P5.
- Target duration per launch (5-10 min proposed) to be set from measurements.
- Whether the roots1 cache (P4) is worth it: depends on measured in-process
  factor-base generation time.
- A test machine: NVIDIA GPU sm_80+, CUDA toolkit, cuda-sieve built
  (`make -C bench`; `fbgen_gpu` for P4).

## Handoff log

Newest first. One entry per session.

### 2026-09-21: P0 + P1 patch applied, builds and runs
- Patch applied to the working tree; yafu compiles and the binary runs.
- Not yet tested against cuda-sieve itself — see "Start here for Phase 2" at
  the top of this file.
- Next: a real cuda-sieve run (c100-c120), then Ctrl-C mid-range, then P2.

### 2026-09-21: P0 + P1 patch drafted
- Produced `cuda-sieve-p0-p1.patch` (7 files, +379/-15, CRLF preserved) and
  `cuda-sieve-p0-p1-tests.tar.gz`. See "Deviations" (9-15) for where it differs
  from the plan.
- Next: apply the patch, do a real yafu build, fix compile issues, then run on a
  GPU machine: c100-c120 end to end, `--check-relations`, msieve filtering,
  Ctrl-C mid-range. Then P2 (salvage from `.ckpt`, exit-code policy) and P3
  (range sizing, multi-GPU).

### 2026-09-21: planning
- Read `nfs_sieving.c`, `nfs.c`, `nfs_impl.h`, `factor.h`, `cmdOptions.{c,h}`,
  `nfs_filemanip.c`, and cuda-sieve's README/RUNBOOK/STATUS/CLI/ckpt/poly code.
- Produced `cuda-sieve-plan.md` and `cuda-sieve-upstream-requests.md`.
- An earlier message in that session described `nfs_sieving.c` before it had
  been opened; it was then read and corrected. Trust the plan, not the early
  wording.
- Next: start P0 (options, fields, sievername helper, `check_for_sievers`),
  then P1.
