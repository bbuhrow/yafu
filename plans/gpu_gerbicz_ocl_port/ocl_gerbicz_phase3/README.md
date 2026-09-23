# OCL_GERBICZ Phase 3 — collision engine

Full OpenCL port of `collision_engine.cu`'s collision-search pipeline:
scatter → per-bucket hash filter → sort/dedup (via Phase 2's radix
sort) → secondary hash → match → emit. Bit-exact against the CPU
reference per decision #7, no sub-group tricks per decision #6.

## Files

| File | Purpose |
|---|---|
| `ocl_collision.h` | Public API. **Read the header comment first** — it explains the base-offset scoping (which of the engine's buffers carry an offset and why), the new `GPU_ARG_LOCAL` arg type, and the overflow-status contract. |
| `ocl_collision_kernels.cl` | All 10 kernels, hand-ported from `collision_engine.cu`. |
| `ocl_collision.c` | Host struct, `ensure_capacity` grow-and-retry loop, `ocl_collision_run()`. |
| `test_ocl_collision.c` | Task 3.5: links directly against the **real Phase 0 harness source** (not a reimplementation) and validates via its own `cmp_run()`. |
| `harness/` | Copy of the Phase 0 harness (`collcase_io.c`, `gen.c`, `cpuref_exact.c`, `cpuref_stats.c`, `cmp.c` + headers) — present in `/mnt/project` already; copied here only so this deliverable builds standalone. |
| `ocl_xface_phase3.h.patch` | Two additions to `ocl_xface.h`: `GPU_ARG_LOCAL` and `GPU_MAX_KERNEL_ARGS` 15→20. Applies cleanly on top of Phase 1's guard patch (verified). |
| `ocl_shared.*`, `ocl_gerbicz_ctx.*`, `ocl_primitives.*`, `ocl_xface_stub.c`, `patched/` | Carried forward unchanged from Phase 1/2. |
| `Makefile` | `make smoke` (gcc), `make CC=clang smoke`, `make mingw-check CL_HEADERS_DIR=<dir>`. |

## Task 3.1 — discrepancies found

1. **`GPU_MAX_KERNEL_ARGS=15` wasn't enough.** Flagged as a *Phase 4*
   concern in Phase 2's STATUS (for the fused trans+scatter kernel) —
   turned out Phase 3's own `ocl_count_and_store_matched_values` already
   needs 16 args. Bumped to 20 (matching `cuda_xface.h`'s existing
   value) rather than picking a new number, so Phase 4 doesn't hit this
   a second time.
2. **CUDA's dynamic shared memory has no OpenCL equivalent in the
   existing `gpu_arg_type_t` enum.** `filter_per_bucket_kernel`'s hash
   table is sized per-launch (`cudaFuncSetAttribute` +
   `MaxDynamicSharedMemorySize`); OpenCL's equivalent is a
   dynamically-sized `__local` kernel *parameter*, set via
   `clSetKernelArg(kernel, idx, size, NULL)` — a fundamentally different
   shape of call than every existing `gpu_arg_type_t` case (which all
   set a *value*, not a *size*). Added `GPU_ARG_LOCAL` for this.
3. **`d_D`/`d_value_counts` don't need separate `_scan`/`_offsets`
   buffers.** CUDA keeps them separate only because
   `cub::DeviceScan::ExclusiveSum` always writes to a distinct output
   buffer; in both cases CUDA never reads the pre-scan buffer again
   once the scanned version exists. Phase 2's `ocl_scan_exclusive_u32`
   scans in place, so this port reuses `d_D` and `d_value_counts`
   themselves as their own post-scan versions — two buffers fewer,
   same semantics.

## A real bug found and fixed during testing

**Cross-queue race.** `ocl_primitives_t` (Phase 2) and this phase's own
`ocl_collision_engine_t` each create their **own separate OpenCL command
queue** via Phase 1's `ocl_thread_init`. Since the two constantly touch
the same buffers back and forth (fill before scatter, reduce-max after
scatter, sort then dedup, scan then scatter_secondary, ...), running
them on two independent queues gives OpenCL no reason to order those
operations relative to each other — a real race, not an efficiency
concern. Symptom: `candidate_count`/`dedup_count` came out plausible
(the scatter→filter→sort→dedup stages all stayed on one side of the
race by luck of enqueue order), but `value_match_count` was **0** in
every single test case. Traced by seeding a sentinel-style debug read
of the secondary hash's bitmask table (`S`) immediately before and
after `count_secondary` ran: it read back as garbage/stale rather than
the freshly-zeroed buffer `ocl_fill_u32` had just written — because
that fill ran on a *different* queue than the kernel that read it.

**Fix applied:** after creating the collision engine's kernel objects
via `ocl_thread_init`, discard the queue it made and reuse
`ocl_primitives_t`'s queue instead, so both halves of the pipeline
share one queue (matching decision #5's "one queue per thread" intent
literally — not one per *kind* of kernel). `ocl_collision_free()`
correctly does not release this borrowed queue.

**Not fixed at the root, flagged for later:** the proper fix is for
`ocl_thread_init()` (Phase 1, `ocl_shared.c`) to accept an *optional*
pre-existing queue instead of always creating one, so callers that want
to share a queue across multiple `ocl_thread_t`s don't need this
create-then-discard workaround. Left as-is this phase to avoid touching
Phase 1's foundation API three phases in; worth doing before Phase 5
has multiple engines/primitive-sets genuinely running concurrently on
separate per-thread queues (where this workaround's assumption — "just
reuse the one other queue that happens to exist" — stops making sense).

## Test results (task 3.5)

Validated against the **real Phase 0 harness** (`gen_generate` +
`cpuref_exact`/`cpuref_stats` for the reference, `cmp_run` for the
comparison — all linked directly, not reimplemented or shelled out to):

| Case | What it exercises | Result |
|---|---|---|
| `tiny_n200` | Harness's own 0.3 edge-case suite (n=200) | found_count=51, exact match |
| `n1_noedge` | n=1, no edge cases | 0 entries, exact match |
| `n5000_density` | n=5000, elevated collision density | found_count=277, exact match |
| `skew500` | 500-key bucket skew — forces `ensure_capacity`'s grow-and-retry loop | found_count=67, exact match |
| `bhash1` | Multiplicative bucket hash (`bucket_hash=1`) | found_count=77, exact match |
| `root8` | `root_bytes=8` (r64) path | found_count=88, exact match |
| `n200000` | n=200,000, several workgroups | found_count=4102, **saturates at 999** — harness's subset-check comparator confirms all 999 stored entries are in the full reference set |
| `smallcap` | `hash_word_cap=32` (small) — forces more filter iterations | found_count=1896, **saturates at 999**, subset check passes |

**8/8 pass**, gcc and clang, zero warnings. The two saturation cases
(`n200000`, `smallcap`) are the ones that actually exercise
`FOUND_ARRAY_SIZE`'s 999-entry cap and the harness's saturation-aware
comparator rule (exact-match below the cap, subset-check at/above it) —
not just the common case. Compile-only mingw-w64 cross-check also
passed with zero warnings on the two new files.

## Task 3.4 — overflow handling

`ocl_collision_status_t` (`OCL_COLLISION_OK` /
`_CANDIDATE_OVERFLOW` / `_VALUE_MATCH_OVERFLOW` / `_INVALID_INPUT` /
`_CL_ERROR`) replaces CUDA's `exit(-1)` on the two overflow paths. Not
wired into `stage1_engine.c`'s `STAGE1_OVERFLOW_SKIP` yet — that's
Phase 5's host-integration job — but the contract exists now so Phase 5
has something concrete to consume rather than inventing one from
scratch.

## Build

```
make smoke                    # gcc, ./patched
make CC=clang smoke
make PROJDIR=/real/path smoke # against the real, both-patches-applied repo tree
make mingw-check CL_HEADERS_DIR=<dir containing only CL/cl.h>
```
