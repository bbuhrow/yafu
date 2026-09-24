# Phase 0 kickoff prompt

**How to use:** start a new conversation in the same Claude Project (so the source files are available under `/mnt/project/`), attach `OCL_GERBICZ_PLAN.md` and `OCL_GERBICZ_STATUS.md`, and paste everything inside the fence below.

**For later phases:** reuse this structure. Keep the "Context", "Working rules", and "End-of-conversation protocol" sections unchanged; replace the phase goal, tasks, source sections, and deliverables with the ones from the plan and the previous phase's handoff notes. The previous conversation should draft that prompt for you.

````
# Context
I'm porting the `gpu_gerbicz` CUDA stage-1 NFS polynomial-selection collision engine (YAFU/msieve) to OpenCL, targeting AMD GPUs that can't run CUDA. The work is split into phases, one per conversation. THIS conversation is **Phase 0 only** (baseline and test harness). No OpenCL or GPU code yet, and don't modify the Project's source files.

Attached: `OCL_GERBICZ_PLAN.md` (full plan, locked decisions, verified facts, line-referenced source map) and `OCL_GERBICZ_STATUS.md` (progress tracker that travels between phases). Read both first. The source files are in `/mnt/project/` — use the line references in the plan's §4 rather than re-reading whole files.

# Working rules
- Keep explanations brief and to the point.
- Deliver work as files, plus one zip of everything. Keep the harness plain C99: no GMP, CUDA, or OpenCL dependencies; no `__int128`; no VLAs; it must build under gcc, clang, and MSVC. Build and test it in your sandbox.
- If a fact in the plan turns out wrong or incomplete, say so and record the correction in STATUS.
- Ask me, in a single message at the start, for anything missing (listed under "Inputs" below), and proceed with sensible defaults meanwhile.

# Phase 0 goal
Produce a portable test harness and reference data so that later phases (the OpenCL collision engine, sieve kernels, and host integration) can be validated bit-for-bit without needing an NVIDIA GPU.

# Tasks
**0.1 Confirm understanding.** Skim the source sections below. List any discrepancy with the plan (short).

**0.2 collcase v1 format.** Define an engine-neutral binary test-case file, little-endian, versioned, with:
- Header: magic, version, `n`, `root_bytes` (4 or 8), `key_bits`, `shift`, `bucket_hash` (0/1), `num_q`, and an optional hash-word cap override.
- Arrays: `keys[n]` (u32 or u64 per `root_bytes`), `values[n]` (u32 = `(q_index << shift) | p`), `q_batch[num_q]` (24-byte `specialq_t`).
- Optional expected-output section: found counter, canonical found entries, and the stats block (`candidate_count`, `dedup_count`, `value_match_count`, `filter_iters_hist[102]`, `bucket_max`).
- Document the byte layout precisely; it will be reused in Phase 3. Include static asserts that `found_t` is 32 B and `specialq_t` is 24 B.

**0.3 Generator.** Deterministic (seeded) synthetic case generator with knobs for `n`, `root_bytes`, `key_bits`, `shift`, `num_q`, and collision density. It must cover:
- signed (negative) keys within ±2^(key_bits−1), stored as 32-bit two's complement or 64-bit signed;
- multi-way collisions, including multiplicity above 8 (exercises the arena fallback);
- equal keys with different `q_index` (must not emit);
- pairs with `gcd(p1,p2) > 1` (must not emit);
- zero keys (skipped);
- bucket skew (exercises bucket growth) with `bucket_hash` 0 and 1;
- p values below and above 65536.

**0.4 CPU exact reference.** Compute the found set at result level, independent of the hash filter (simple sort/group). Semantics follow the CUDA pipeline: skip `k == 0`; for each key group, pair entries whose `q_index` is equal and `gcd(p1,p2) == 1`; the emitted entry is `{p1, p2, q = q_batch[qi].p, qroot = q_batch[qi].root, offset = signed key}`; root sign rule for 4-byte roots is `(int64)(int)key`. Report the total attempted-store count (the value the CUDA found-array counter would reach) and the entries.

**0.5 CPU stats model.** A sequential, exact mirror of the filter pipeline (scatter to `NUM_BUCKETS` buckets via `compute_bucket`, per-bucket filter, dedup, secondary hash, matching) that reproduces `candidate_count`, `dedup_count`, `value_match_count`, `filter_iters_hist[102]`, and `bucket_max`, given a configurable hash-word cap (default: derived from a 64 KB local-memory limit, i.e. 4096 words). The filter is order-independent, so an exact match is possible. Cross-check that the model's found set equals the 0.4 reference. If this task slips, record that in STATUS: Phase 3b's exit criteria depend on it.

**0.6 Comparator.** A `cmp` mode that canonicalizes results (order `p1 < p2` within an entry; sort by `(q, qroot, offset, p1, p2)`), compares two result files, and reports differences. Saturation rule: compare full sets only when the counter is below 999; otherwise compare the counter and check that stored entries are a subset of the reference.

**0.7 Named suite.** A script that deterministically generates a small named suite, e.g. `tiny_*` for unit tests, `med_*` (about 1–4M elements), and one optional large (about 30M) case; cover root_bytes 4 and 8 and the edge cases from 0.3. Include the reference outputs.

**0.8 End-to-end recipe.** Write down exact steps to produce `test_cpu_gerbicz.hits` for chosen `test_ad` cells (`nfs_args`: `test_ad= test_pmin= test_pmax= test_qmin= test_qmax=`), how to produce the same file with the OpenCL engine later (`test_ocl_gerbicz.hits`), and the comparison rule (sort lines, diff; hits are `a_d p m`, post-filter, single thread). Propose a cell list with me: small input (<480 bits), a larger input, pp32 (p_max < 65536) and pp64 (p_max ≥ 65536).

**0.9 (Optional, only if I have an NVIDIA machine)** Specify a small patch (do not apply it to the Project) for `collision_engine_run` that dumps a collcase v1 file with the CUDA engine's actual outputs and stats, for use as golden data.

**0.10 Wrap-up.** Update STATUS (phase table, interfaces and formats, artifacts with filenames, open issues, handoff notes) and draft the Phase 1 prompt using this file's structure.

# Source sections to read (from `/mnt/project/`)
- `collision_engine.cu`
  - 30–37 constants
  - 53–76 `compute_ilog2`/`compute_capped_ilog2`
  - 94–153 `store_hit_collision` and `scatter_roots_kernel`
  - 155–336 `filter_per_bucket_kernel`
  - 338–407 dedup, secondary, `find_candidate_slot`
  - 409–598 match and emit kernels
  - 685–758 `ensure_capacity` (bucket sizing)
  - 807–814 `host_ilog2`
  - 866–1162 `collision_engine_run` (pipeline, hash cap at 928–961, overflow at 1005–1013 and 1108–1113, arena decision at 1073–1077)
- `collision_bucket.h` (whole file)
- `collision_engine.h` 14–60 (`collision_data_t`)
- `stage1_core.h` 26–42 (`FOUND_ARRAY_SIZE`, `found_t`, `specialq_t`)
- `stage1_sieve_gpu.c`
  - 611–662 `check_found_array`
  - 780–831 collision call and field setup
  - 1064–1106 batching, `key_bits`, `num_aprog_vals`, `shift`
- `stage1.c` 118–200 (`handle_collision`, test-mode dump), 522–551 (test-mode args)

# Inputs (ask me in one message)
1. Upload `ocl_xface.h` to the Project (needed from Phase 1).
2. A representative YAFU invocation and input N (and `nfs_args`) for the end-to-end `test_ad` cells.
3. Do I have an NVIDIA machine available for optional CUDA dumps (task 0.9)?
4. Where in my repo should the harness live (path convention)?

# End-of-conversation protocol
Deliver: the harness sources, the suite generator script, a README (build, run, format spec), a zip of everything, the updated `OCL_GERBICZ_STATUS.md`, and the Phase 1 prompt. Then give me a 3–5 line summary of what was done, what was deferred, and anything I need to do before Phase 1.
````
