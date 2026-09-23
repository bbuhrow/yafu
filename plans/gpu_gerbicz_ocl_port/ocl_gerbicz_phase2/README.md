# OCL_GERBICZ Phase 2 — primitives

Implements the four building blocks Phase 3's collision engine needs in
place of CUB (`cub::DeviceReduce::Max`, `cub::DeviceRadixSort::SortKeys`,
`cub::DeviceScan::ExclusiveSum` x2 — `collision_engine.cu` 741–757):
buffer fill, reduce-max, exclusive scan, and an LSD radix sort over
64-bit keys. Built on Phase 1's `ocl_shared.h`/`ocl_gerbicz_ctx.h`.

## Files

| File | Purpose |
|---|---|
| `ocl_primitives.h` | Public API. **Read the header comment first** — it explains why this phase implements only ONE radix sort (over `ulong`), not separate 32-/64-bit variants, and the two-level scan's documented ~16.7M-element capacity limit. |
| `ocl_primitives_kernels.cl` | The eight kernels (2 fill, 1 reduce-max, 3 scan, 2 radix). Chunk-size constants arrive as `-D` flags from `ocl_primitives.c` so host and device can't drift apart. |
| `ocl_primitives.c` | Host wrappers: program build (via Phase 1's `ocl_build_program_cached`), scratch-buffer management, the reduce-max/scan/radix-sort host-side loops. |
| `test_ocl_primitives.c` | Task 2.6/2.7: CPU references, correctness tests across edge sizes through CANDIDATE_CAP scale, and rough benchmarks. |
| `ocl_shared.*`, `ocl_gerbicz_ctx.*`, `ocl_xface_stub.c`, `patched/` | Carried forward from Phase 1 unchanged (per your answer to question 3 — alongside for now, you'll distribute later). |
| `Makefile` | `make smoke` (gcc), `make CC=clang smoke`, `make mingw-check CL_HEADERS_DIR=<dir>`. |

## Task 2.1 — discrepancies found

1. **The plan's "~4M" scan size and `MAX_DSIZE` are two different things.**
   `MAX_DSIZE` (`collision_engine.cu:36`) is `(1u<<20)+64` ≈ **1.05M**, not
   ~4M. The actual ~4M `cub::DeviceScan::ExclusiveSum` call is a
   *different* one, three lines later (`collision_engine.cu:754-756`,
   over `d_value_counts`/`d_value_offsets`, size `VALUE_MATCH_CAP+1` =
   `CANDIDATE_CAP+1` = **4,194,305**). So the plan's phrasing wasn't
   wrong — it just needed disambiguating between two call sites. Tested
   both sizes explicitly (see test list below).
2. **Only one radix sort is needed, not two.** `scatter_roots_kernel`
   (`collision_engine.cu:122-126`) assigns a 4-byte root to a `uint64 k`
   — an implicit *zero*-extension (sign-extension only happens later, on
   emit). `d_candidate_keys`/`d_sorted_keys` are declared `uint64*`
   unconditionally regardless of `root_bytes`. So there is only ever one
   storage width at the point sorting happens; "sort the full 32- or
   64-bit key, not `key_bits`" (the plan's Section 3 fact) is about never
   truncating to `key_bits`, not about needing two sort implementations.
   Implementing only the 64-bit sort avoids duplicating a full working
   sort for a case the codebase never actually exercises. Also follows
   from this: the sort can be a plain *unsigned* LSD radix sort (no
   sign-bit-flip trick) — matches `cub::DeviceRadixSort::SortKeys`'s own
   default behavior on a plain `uint64_t*` buffer.
3. **`GPU_MAX_KERNEL_ARGS` wasn't bumped in the OpenCL header when the
   CUDA one was.** `cuda_xface.h` (now attached) shows it was bumped
   15→20 on 2026-05-25 for the fused trans+scatter kernel's 17 args, but
   `ocl_xface.h`'s copy is still `15`. Not a problem for Phase 2 (largest
   kernel here uses 9), but **Phase 4's fused trans+scatter port will
   need this bumped** — flagging now so it isn't a surprise then.

## A real bug found and fixed during testing (worth calling out)

`ocl_scan_local_u32`'s block-sums write (`block_sums[wg] = acc;`) was
present in the design comment above the kernel but **missing from the
actual kernel code** — a plain omission, not a subtlety. It went
undetected until the test harness exercised more than one block-scan
workgroup (`n > 4096`): with exactly one workgroup, the missing write
doesn't matter (the single block's "block sum" is always scanned down
to `0` anyway), so the `n=4096` test passed while masking the bug, and
`n=4097` caught it immediately with garbage in every element past the
first block. Found by bisecting with a sentinel-fill of the scratch
buffer before the kernel ran (confirmed the kernel wasn't touching
`block_sums` at all, not just computing it wrong) — see the git history
of `ocl_primitives_kernels.cl` if you want the exact diff. Fixed, and
now covered by `n=4097`/`n=8193`/`n=100000` and the two large sizes in
the test list, all of which exercise the multi-workgroup path.

A second bug, also found via testing: the radix scatter kernel's
per-workgroup bin cursor (`cursor[PRIM_RADIX_BINS]`, 256 `uint`) was
originally a **private** (per-work-item) array. Even though only
`lid==0` ever uses it, PoCL's CPU backend appears to size private memory
per work-item across the whole work-group regardless of whether a given
lane's branch reaches the declaration, and 256 work-items × 1KB blew the
per-work-group stack, segfaulting only at larger scale
(`num_wg` ≥ roughly 25 in this sandbox). Moved to `__local` (one
workgroup-shared allocation, still only touched by `lid==0`) and the
crash disappeared. Worth remembering for any future kernel that declares
a moderately large array and returns early on `lid!=0`: prefer `__local`
over a private array unless every lane genuinely needs its own copy.

## Test results (task 2.6)

All sizes below tested on both fill/reduce-max/scan and (except fill)
radix sort where applicable: `0, 1, 17, 255, 256, 4096, 4097, 8193,
100000`, plus the real sizes actually used in `collision_engine.cu`:
`16384` (`NUM_BUCKETS`, reduce-max), `1048640` (`MAX_DSIZE`, scan),
`4194305` (`VALUE_MATCH_CAP+1`, scan), `4194304` (`CANDIDATE_CAP`, radix
sort — both uniform-random and an adversarial distribution with heavy
duplication and the zero-extended-negative-root bit pattern from the
plan's Section 3 fact).

**41/41 passed**, gcc and clang, zero warnings under
`-Wall -Wextra -pedantic`. Compile-only mingw-w64 cross-check also
passed with zero warnings (see Phase 1's README for what that
does/doesn't prove — still not a substitute for real MSVC).

## Benchmarks (task 2.7)

You confirmed a real AMD card is reachable — but reachable to **you**,
not to this sandbox: this container has no `/dev/dri`, no ROCm, and no
route to your machine, so the numbers below are the same PoCL
CPU-emulated OpenCL device as Phase 1, not your RX 6700 XT. Treat them
as "the pipeline runs and produces a number", not as anything
predictive of real GPU throughput:

```
scan_exclusive_u32   n=4,194,305   0.015 s   (~288 M elems/s, CPU-emulated)
radix_sort_u64       n=4,194,304   0.46  s   (~9 M keys/s,   CPU-emulated)
```

**To get a real number**, copy this deliverable onto the machine with
the RX 6700 XT, `make CL_HEADERS_DIR=<your AMD OpenCL SDK's include dir
if needed> smoke` (or just `make smoke` if your system's OpenCL headers
are already on the default include path), and the benchmark lines will
print the same way. Send me that output and I'll fold it into STATUS —
I can't run it there myself.

Known, deliberate non-optimizations that a real benchmark will show up
(all flagged "Phase 7 target" in the source comments, per decision #7
"bit-exact first, then optimize"):
- The scan's per-block combine step (`thread_sum[256]` scan) is done by
  a single thread doing 256 serial adds, not a parallel Hillis-Steele
  scan.
- The radix scatter kernel is single-thread-per-workgroup (one lane
  walks its whole 4096-element chunk sequentially) rather than using a
  proper parallel local prefix sum.
- Radix width is 8 bits/pass (8 passes for a 64-bit key) with no
  attempt made to tune this against the target hardware yet (task 2.5
  said no preference given, so picked the common default).

## Build

```
make smoke                    # gcc, ./patched (Phase 1's guard-fixed copy)
make CC=clang smoke
make PROJDIR=/real/path smoke # against the real, ocl_xface.h.patch-applied repo tree
make mingw-check CL_HEADERS_DIR=<dir containing only CL/cl.h>
```
