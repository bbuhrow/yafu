# OCL_GERBICZ Phase 1 — OpenCL foundation

Delivers the OpenCL context/program/queue/kernel infrastructure that
Phase 3+ kernels will run inside. No collision-engine kernels are ported
here (see OCL_GERBICZ_PLAN.md).

## Files

| File | Purpose |
|---|---|
| `ocl_shared.h/.c` | Generic device/context/program-cache/thread infrastructure. Not gerbicz-specific — usable by `gpu_cofactorization_cl.c` too. |
| `ocl_gerbicz_ctx.h/.c` | Gerbicz-specific layer: two-program split, sub-group build flag, `-D` flags from `collision_bucket.h`. **Read the big comment block at the top of `ocl_gerbicz_ctx.h` first** — it's the base-offset argument convention and sub-group gating design, the two things every later kernel signature depends on. |
| `ocl_xface.h.patch` | One-line-guard fix for the real `ocl_xface.h` (see "Discrepancies found" below). Apply to the real project tree with `patch -p1 -d <path> < ocl_xface.h.patch`. |
| `ocl_xface_stub.c` | **Sandbox-only.** Reference implementation of `gpu_init`/`gpu_launch_init`/`gpu_launch_set`, which are declared in `ocl_xface.h` and already called by `gpu_cofactorization_cl.c`, but whose implementation wasn't among the files handed to this phase. Needed only so this deliverable's own smoke test can build and run; see "Open issues" before relying on it for anything real. |
| `ocl_gerbicz_smoke.cl` / `test_ocl_foundation.c` | Task 1.7 smoke test: trivial kernel exercising the base-offset convention and the `HAVE_SUBGROUPS` gate, run through the full build→cache→launch→verify pipeline twice (compile, then cache-hit). |
| `patched/` | A copy of `ocl_xface.h` (patched) and `collision_bucket.h` (unmodified), so `make smoke` works out of the box without touching the real project tree. Point `PROJDIR` at the real tree instead once the patch is applied there. |
| `Makefile` | `make smoke` (gcc), `make CC=clang smoke`, `make mingw-check CL_HEADERS_DIR=<dir-with-only-CL/cl.h>`. |

## Discrepancies found against the plan (task 1.1)

1. **`ocl_xface.h`'s feature guard would silently empty the header for a
   gerbicz-only build.** It's `#ifdef HAVE_OCL_BATCH_FACTOR` (the
   cofactorization flag). Locked decision #8 proposes a separate
   `HAVE_OCL_POLY` guard for gerbicz; building with only that defined
   would make every declaration in `ocl_xface.h` disappear (the include
   guard becomes a no-op), not fail loudly. Fixed by
   `ocl_xface.h.patch` (relaxes the guard to `defined(HAVE_OCL_BATCH_FACTOR)
   || defined(HAVE_OCL_POLY)`). No locked decision needed to change —
   decision #3 already says the two OpenCL builds aren't required to be
   mutually exclusive of each other, only of CUDA.

2. **`gpu_info_t` has no `MAX_MEM_ALLOC_SIZE` field**, though the plan's
   Phase 1 bullet list asks for it. `shared_mem_size` and
   `max_threads_per_block` cover local-mem-size and max-work-group-size,
   but there's no equivalent for `CL_DEVICE_MAX_MEM_ALLOC_SIZE`. Not
   worth a locked-decision change or a patch to a struct that other code
   already depends on — `ocl_device_t` (this phase) just queries it
   directly from `cl_device_id` and stores it alongside the `gpu_info_t*`
   it wraps.

3. **Device version is not proof a given `-cl-std=` builds.** Found
   empirically in this sandbox: PoCL reports `CL_DEVICE_VERSION` =
   "OpenCL 3.0" but rejects `-cl-std=CL2.0` outright — its
   `CL_DEVICE_OPENCL_C_ALL_VERSIONS` list is `1.0, 1.1, 1.2, 3.0`, with
   no 2.0 entry at all. `ocl_device_query()` now queries that list (with
   a fallback to the older single-string `CL_DEVICE_OPENCL_C_VERSION`
   query) and stores the actual best usable `-cl-std=` value in
   `ocl_device_t::cl_c_std`, rather than trusting `compute_version_major/minor`
   for this purpose. This doesn't change the locked-decision #1 *floor*
   (still refuse < 2.0, via `ocl_check_min_version` on the device
   version) — it only changes how the build-option string is picked once
   that floor is cleared. Worth remembering for the real AMD target too:
   don't hardcode `"-cl-std=CL2.0"` anywhere in Phase 3+.

4. `gpu_init`, `gpu_launch_init`, `gpu_launch_set` are declared in
   `ocl_xface.h` and already called from `gpu_cofactorization_cl.c`, but
   no implementation file for them was among this phase's inputs. Not
   necessarily a bug (it may simply live in a file not handed to this
   phase) — flagged so it gets reconciled before Phase 3 builds on it.
   See `ocl_xface_stub.c` and "Open issues" below.

Everything else in the plan/STATUS matched the source as read.

## Locked-decision implementation notes

- **Decision #2 (cache):** fixed the Phase 0-flagged weakness in
  `gpu_cofactorization_cl.c`'s own cache (device-name-only key) —
  `ocl_build_program_cached()` keys on device name + driver version +
  build-options string + an FNV-1a hash of the concatenated source
  text. Verified in this sandbox: editing the kernel source correctly
  forces a recompile rather than silently reusing a stale binary (see
  "Smoke test results" below).
- **Decision #4 (base-offset args):** the full convention (element
  offsets, always `cl_ulong`/`ulong`, argument ordering, when an offset
  is/isn't needed) is written up in `ocl_gerbicz_ctx.h`'s header
  comment. This needed to be pinned down now since every Phase 2/3
  kernel signature depends on it.
- **Decision #5 (context/program/queue split):** `ocl_device_t` holds
  one context (+ per-purpose `cl_program`s in `ocl_gerbicz_device_t`);
  `ocl_thread_t` holds one queue + one `cl_kernel` per entry point,
  matching `gpu_cofactorization_cl.c`'s existing per-thread pattern and
  `collision_engine`'s own init/free lifecycle shape (mirrored, not
  reused directly, since Phase 3 owns the actual buffers).
- **Decision #1/#6 (sub-groups):** build-time gate, not runtime branch —
  see the rationale in `ocl_gerbicz_ctx.h`. `ocl_gerbicz_build_opts()`
  sets `-DHAVE_SUBGROUPS=0/1` from a real extension-string query
  (`cl_khr_subgroups`/`cl_intel_subgroups`), plus bakes in
  `LOG2_NUM_BUCKETS`/`BUCKET_HASH_MIX` from `collision_bucket.h` as
  further `-D` flags (per the plan's Phase 1 bullet) rather than letting
  `.cl` sources redefine them and risk drifting from the header that
  explicitly warns against exactly that.

## Smoke test results (task 1.7)

Ran in this sandbox (no AMD hardware available here either — see "no
NVIDIA machine" precedent from Phase 0). Installed a CPU OpenCL runtime
(PoCL 5.0) plus the ICD loader and Khronos headers, since none were
present in the base image; also installed clang-16 (only gcc-13 was
present) and mingw-w64 (for the compile-only Windows check).

- **gcc 13 and clang 18(-equivalent, clang-16 here)**: `make smoke` /
  `make CC=clang smoke` both PASS, zero warnings under
  `-Wall -Wextra -pedantic`.
- **Pipeline exercised end-to-end twice per run**: first launch compiles
  from source and writes the binary cache; second launch hits the cache
  (prints `loaded cached binary`). Both verified correct against 1024
  elements.
- **Base-offset convention actually exercised**, not just declared: all
  three buffers (`a`, `b`, `out`) are allocated larger than needed and
  read/written from non-zero offsets (7, 13, 3 elements) baked into the
  kernel launch, so a convention bug would show up as wrong output.
- **Cache invalidation verified**: deliberately edited the kernel's
  arithmetic and reran — the changed source was correctly recompiled
  (not silently served from the stale cache) and the resulting mismatch
  was correctly detected and reported (`FAIL`, `got 1001 want 1`),
  confirming the source-hash-keyed cache actually protects against the
  Phase 0-flagged bug. This is also incidentally why the sub-group path
  is real: this device reports `cl_khr_subgroups` support, so
  `HAVE_SUBGROUPS=1` was actually compiled and run, not skipped.
- **MSVC**: not attempted. No `cl.exe` toolchain is reachable from this
  sandbox (Linux container; the allowed network egress list has no
  Microsoft domains). As a partial substitute, ran a **compile-only**
  cross-check under mingw-w64 (a GNU toolchain targeting Windows,
  `make mingw-check`) — this passed with zero warnings, which rules out
  gross portability breaks (missing headers, GNU-only extensions,
  reliance on POSIX-only APIs) but is **not** a substitute for MSVC
  itself: it can't catch MSVC-specific issues (older MSVC's C89-only
  declaration placement rules, `snprintf`/`sscanf` conformance
  differences, `_CRT_SECURE_NO_WARNINGS`-class warnings). A real MSVC
  build is still owed the first time this code is touched from Windows.

## Open issues carried forward

- `ocl_xface_stub.c` is a stand-in. If the real `ocl_xface.c` already
  exists elsewhere in the repo, replace the stub and re-run the smoke
  test against the real implementation before Phase 3 relies on any of
  this — in particular check that its `gpu_launch_set` handles the same
  five `gpu_arg_type_t` cases the same way, and that a real
  AMD `warp_size`/preferred-multiple isn't hardcoded to 32 anywhere (the
  stub's `gpu_init` hardcodes `gi->warp_size = 32` as a placeholder,
  since getting the real value needs a built kernel — `ocl_device_t`'s
  own `pref_wg_multiple` is the value to actually use, not this field).
- `clGetErrorString`'s `default:` case uses a `static char buf[32]` —
  not thread-safe. Fine for Phase 1 (single-threaded smoke test); flag
  before any multi-threaded phase (5+) calls `OCL_TRY` concurrently from
  two threads on an unrecognized error code.
- MSVC build still unverified for real (see above).
