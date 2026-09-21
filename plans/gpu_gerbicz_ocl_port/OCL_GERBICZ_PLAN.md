# OpenCL port of the `gpu_gerbicz` stage 1 engine — plan

Target: AMD GPUs that cannot run CUDA (dev card: RX 6700 XT, RDNA2, OpenCL 2.0).
Scope: Gerbicz collision engine only. **No sort-engine port.**
Work is split into phases, one phase per conversation. Progress lives in `OCL_GERBICZ_STATUS.md`.

---

## 1. Scope

| Piece | CUDA source | Port? |
|---|---|---|
| Collision engine (scatter, filter, sort/dedup, secondary hash, match, emit) | `collision_engine.cu/.h`, `collision_bucket.h` | Yes |
| Sieve kernels `sieve_kernel_trans_pp32_r32 / pp32_r64 / pp64_r64` | `stage1_core.cu` 22–613 | Yes |
| `sieve_kernel_final_32/64`, `store_hit` | `stage1_core.cu` 615–711 | **No** (sort-engine path only) |
| Host driver | `stage1_sieve_gpu.c` | Yes, as new `stage1_sieve_ocl.c` |
| Intrinsics | `cuda_intrinsics.h` | Already ported: `opencl_intrinsics.cl` |
| Engine registry | `stage1_engine.c/.h` | Add one row |

## 2. Locked decisions

1. Minimum OpenCL 2.0 (`-cl-std=CL2.0`). Work-group functions are allowed. Sub-groups are an optional fast path, gated on the extension.
2. Static link, `.cl` files loaded at runtime. No DSO (the CUDA DSO exists only for nvcc/CUB).
3. New symbols are `ocl_`-prefixed. CUDA and OpenCL builds are mutually exclusive for now.
4. Pointer offsets: kernels get base-offset arguments (`p_out += off` inside the kernel). No sub-buffers, no SVM.
5. One context and program per device (shared). One queue and its own kernel objects per host thread (`clSetKernelArg` is not thread-safe).
6. Baseline collision kernels use no sub-groups: local and global atomics only (the existing `SCATTER_DIRECT_ATOMIC` path is the model).
7. Bit-exact first: port verbatim, validate, then optimize. Behavior changes (bug fixes) are made on both backends or logged in STATUS.
8. Proposed registry entry: id `STAGE1_ENGINE_OCL_GERBICZ`, name `ocl_gerbicz`, guard `HAVE_OCL_POLY` (guard name to be confirmed in Phase 5).

## 3. Facts established from the source

- **Keys are signed.** The trans kernels do `if (newroot > pp/2) newroot -= pp`. Roots are stored as 32-bit two's complement (r32) or 64-bit signed (r64). The engine zero-extends 4-byte roots to 64 bits (`collision_engine.cu` 122–126) and sign-extends on emit (541–542). So the radix sort must sort the full 32- or 64-bit key, not `key_bits`. Sorting only `key_bits` bits is a later optimization, valid only if masking to `key_bits` is injective on the real range.
- **Key 0 means "empty".** Scatter and match kernels skip `k == 0`. The trans kernels write 0 for non-coprime (p, q) pairs.
- **The packed value is `(q_index << shift) | p`.** `q` here is the index into `q_batch`, not the q value. Emit pairs entries only when the indices are equal and `gcd(p1,p2) == 1`. The q value and root come from `q_batch[q_index]`.
- **Found array.** `found_array[0].p1` is an atomic counter of all attempted stores. Entries are stored only while `index < FOUND_ARRAY_SIZE-1` (999 usable). Sets are comparable only when the counter is below that.
- **Overflow handling.** The registry has `STAGE1_OVERFLOW_SKIP` but nothing consumes it. The CUDA engine calls `exit(-1)` on candidate overflow (`collision_engine.cu` 1005–1013) and value-match overflow (1108–1113).
- **Batch cap.** `stage1_sieve_gpu.c` 1051–1057 hard-caps the special-q batch at 16384 for the collision engine.
- **Hash-table cap.** The number of shared-memory hash words comes from the opt-in shared-memory limit (`collision_engine.cu` 941–961). With 64 KB of OpenCL local memory that caps at 4096 words. Watch `hash_cap_count`.
- **Filter is order-independent.** The per-bucket filter result (candidate set and count, stop-iteration histogram) depends only on bucket contents and parameters. A sequential CPU model can match it exactly.
- **Quirk in `pp64_r64`.** `stage1_core.cu` line 502: `uint32 curr_qq_prod = roots_out[...]` truncates a 64-bit value (line 517 uses `uint64`). Port verbatim first.
- **Trans kernel scratch.** `roots_out` (r32, pp64_r64) or `p_out` (pp32_r64) is used as scratch for qq_prod, then overwritten. For `num_aprog_vals > 1` the host pre-clears the roots array (`stage1_sieve_gpu.c` 689–700).
- **Test mode.** `nfs_args`: `test_ad= test_pmin= test_pmax= test_qmin= test_qmax=`. It forces 1 thread and writes `test_<engine>.hits`, one line per hit: `a_d p m`, post-filter (`stage1.c` 195–198, 528–551, 674–690, 915–923).
- **Precedent** (`gpu_cofactorization_cl.c`):
  - The cache file key is the device name only (1344–1357).
  - Build options are `-cl-std=CL2.0 -cl-mad-enable` (1340).
  - One context and program per thread (1295+), in-order queue without profiling (1501), wall-clock timing.
  - `ocl_xface.h` (not in the Project yet) supplies `gpu_info_t`, `gpu_launch_t`, and `OCL_TRY`.
- **Intrinsics signature differences** (`opencl_intrinsics.cl`): `modinv64(a, p, ulong *likely_gcd)` (line 500) and `montmul64(a, b, n, uint w)` (809/838) differ from the CUDA signatures.

## 4. Reference map

**`collision_engine.cu`** (1164 lines)

| Lines | Content |
|---|---|
| 30–37 | Constants: MAX_FILTER_ITERS 20, BLOCK_THREADS 128, CANDIDATE_CAP 2^22, MATCH_ARENA_WIDTH 8, MAX_C_ILOG2 20 |
| 53–76 | `compute_ilog2`, `compute_capped_ilog2` |
| 94–109 | `store_hit_collision` |
| 111–153 | `scatter_roots_kernel` (131–146: direct-atomic vs `__match_any_sync` path) |
| 155–336 | `filter_per_bucket_kernel` (198–206 first pass; 212–317 iteration loop; 319–335 emit) |
| 338–385 | dedup, count_secondary, scatter_secondary |
| 387–407 | `find_candidate_slot` |
| 409–524 | The three match kernels |
| 526–598 | `emit_found_kernel`, `emit_found_arena_kernel` |
| 600–805 | `struct collision_engine`; `ensure_capacity` at 685–758 |
| 866–1162 | `collision_engine_run`: 886–926 scatter and grow loop; 928–961 hash cap; 963–1013 filter and overflow; 1015–1037 sort and dedup; 1039–1064 secondary hash; 1066–1112 match; 1114–1161 emit and fallback |

**Other files**

- `collision_bucket.h` 1–40 (whole file): constants and `compute_bucket`.
- `collision_engine.h` 14–60: `collision_data_t` (uses `CUdeviceptr`/`CUstream`; must become `cl_mem` and `cl_command_queue`).
- `stage1_core.h` 26–42: `FOUND_ARRAY_SIZE`, `found_t` (32 B), `specialq_t` (24 B).
- `stage1_core.cu`: `pp32_r32` 22–215, `pp32_r64` 218–415, `pp64_r64` 418–613.
- `stage1_sieve_gpu.c`:
  - 126–359 `p_soa` arrays
  - 363–458 `specialq` arrays
  - 460–551 thread and device structs
  - 559–598 emergency cleanup
  - 611–662 `check_found_array`
  - 668–886 `handle_special_q_batch` (trans launches 709–778, collision branch 780–831)
  - 908–1320 `sieve_specialq` (sizing 985–1001, batch cap 1045–1057, batching and key_bits 1064–1106)
  - 1324–1414 `stage1_specialq_gpu`
  - 1603–1739 thread init/free
  - 1742–1847 device init
- `stage1_engine.h/.c`: registry (CUDA rows at `.c` 68–86).
- `stage1.c`: `handle_collision` 118+, engine selection 522, test mode 528–551, `is_gpu` at 208 and 757.
- `opencl_intrinsics.cl`: `gcd32` 444, `modinv32` 469, `modinv64` 500, `montmul32` 780, `montmul64` 809, `montmul64_r` 949.
- `gpu_cofactorization_cl.c`: `gpu_ctx_init` 1295–1513, cache 1338–1481, `gpu_ctx_free` 1522.

## 5. Phases

### Phase 0 — Baseline and harness
- Engine-neutral test-case file format ("collcase").
- Generator: synthetic cases with signed keys, planted collisions, multi-way keys, cross-q duplicates, gcd>1 pairs, zero keys, bucket skew, root_bytes 4 and 8.
- CPU exact reference (found set and counter), independent of the hash filter.
- CPU stats model that mirrors the filter pipeline exactly: `candidate_count`, `dedup_count`, `value_match_count`, `filter_iters_hist[102]`, `bucket_max`, with an overridable hash-word cap.
- Comparator with the canonicalization and saturation rules.
- Named test-case suite.
- End-to-end reference recipe: `test_ad` cells run on `cpu_gerbicz`.
- Optional: a CUDA dump-hook spec, if an NVIDIA machine exists.

**Exit:** harness builds and self-checks (reference vs stats model agree); the suite is generated; the end-to-end recipe and cell list are written down.

### Phase 1 — OpenCL foundation
- Read `ocl_xface.h` (required input).
- Extract the shared code from `gpu_cofactorization_cl.c`: device selection, program build, cache, `OCL_TRY`. Cofactorization must keep building on it.
- Fix the cache key: source hash + driver version + build options + device.
- One program per device, per-thread queues and kernels.
- Query and expose `LOCAL_MEM_SIZE`, `MAX_MEM_ALLOC_SIZE`, `MAX_WORK_GROUP_SIZE`, and the preferred work-group multiple. Do not assume a wave size of 32.
- Generate `-D` options from `collision_bucket.h`.
- Separate programs for the sieve kernels and the collision kernels.

**Exit:** trivial kernels run per thread on the target card; cofactorization unaffected.

### Phase 2 — Primitives (replacing CUB)
- Exclusive scan (up to about 4M), reduce-max (16384), buffer fill, LSD radix sort of 32- or 64-bit keys (full key width; see Section 3).
- Block-local scan plus block-sums pass; no decoupled look-back.
- Tests against CPU references, plus benchmarks on the target card.

**Exit:** randomized and edge-size tests pass.

### Phase 3 — Collision engine
- **3a:** scatter, reduce-max, host grow-and-retry loop.
- **3b:** per-bucket filter.
  - 128-thread groups with `reqd_work_group_size`.
  - Survivor compaction via a local atomic counter (output order does not matter).
  - Hash tables sized from `LOCAL_MEM_SIZE`.
- **3c:** sort, dedup, secondary hash, count/store/scatter, emit.
- **3d:** fallbacks (arena overflow, `VALUE_MATCH_CAP`, found-array saturation).
- Replace `exit(-1)` with a return code that the driver treats as a skip.
- Explicit `found_t`/`specialq_t` layouts, with host-side size asserts.
- Standalone harness runs against the collcase suite.

**Exit:** found sets and counts match the CPU reference; the stats match the model when the hash cap is aligned.

### Phase 4 — Sieve kernels
- Port the three trans kernels with base-offset arguments, using `opencl_intrinsics.cl`.
- Reconcile the `modinv64`/`montmul64` signature differences.
- Verbatim first (line 502 quirk).
- Validate against a small host-side reference on sample (p, q) pairs, then end-to-end in Phase 6.

**Exit:** roots match the reference for all three variants, including `num_aprog_vals > 1`.

### Phase 5 — Host integration
- New `stage1_sieve_ocl.c`, registry row, and the SKIP policy wired into `stage1.c`.
- Port `p_soa`, `specialq`, `check_found_array`, and the batching loop.
- Launch geometry: global `(blocks_x*size_x, blocks_y)`, local `(size_x, 1)`. `size_x` steps by the preferred multiple, not `warp_size`. Retune the `total_blocks` heuristic (`stage1_sieve_gpu.c` 729–753).
- Profiling or wall-clock timing, keeping the 60 s hibernation guard.
- Cap allocations by `MAX_MEM_ALLOC_SIZE`; redo the atexit cleanup.

**Exit:** a full `test_ad` run completes.

### Phase 6 — Validation
- Diff sorted `test_ad` hit sets against `cpu_gerbicz` (and CUDA where available).
- Matrix: degrees 4–6, pp32/pp64, r32/r64, `num_aprog_vals > 1`, small inputs (<480 bits), threads 1–4, `which_gpu`.
- Log saturation and skipped cells; soak run.

**Exit:** hit sets match, apart from documented saturation and skips.

### Phase 7 — Performance and hardening
- Filter occupancy (local memory per group); wave32 vs wave64.
- Global-atomic contention on the 16384 bucket counters.
- Optional sub-group scatter and compaction.
- Cutting the roughly 5 host round-trips per batch.
- Coalescing and the 64-bit `%` cost in the trans kernels.
- Embedded kernel sources, build flags, docs.

## 6. Working agreements

- One phase per conversation. Attach this plan and the STATUS file at the start.
- At the end of each phase: update STATUS (phase table, decisions or facts changed, artifacts with filenames, interfaces, open issues) and draft the next phase's prompt in the same structure as `OCL_GERBICZ_PHASE0_PROMPT.md`.
- Deliverables are files (zip the bundle). Project sources are not modified in place; changes are delivered as new files or patches.
- Keep explanations brief.
