# GPU `gpu_gerbicz` Collision Engine — Architecture Reference

Companion to `poly_architecture.md`, scoped to stage-1's GPU collision engine
(`stage1_engine=gpu_gerbicz`, engine id #4): the host-side batching in
`stage1_sieve_gpu.c` and the device kernel pipeline in `collision_engine.cu`.
Written after diagnosing a `candidate overflow` crash; includes that failure
mode and its fix as a first-class gotcha (§7).

---

## 1. What it does, in one breath

Same collision search as the other stage-1 engines (find `(p1,p2,root)`
triples that collide on a special-q's arithmetic progressions), but instead
of a global sort it **hash-buckets** all generated roots, then iteratively
narrows each bucket toward its true duplicates using progressively wider
slices of the key's bits, and finally cross-references survivors against the
per-thread (p,q) value stream to emit `found_t` hits.

Two source files own this:

| file                  | role                                                        |
|-----------------------|-------------------------------------------------------------|
| `stage1_sieve_gpu.c`  | host driver — builds arithmetic progressions, batches special-q, derives `key_bits`, calls the engine |
| `collision_engine.cu` | device kernels + `collision_engine` C++ class (`ensure_capacity`, buffer lifetime) — the actual bucket/filter/dedup/match/emit pipeline |

---

## 2. Host side: batching and `key_bits` (`stage1_sieve_gpu.c`)

### 2.1 Arithmetic progressions and special-q

`sieve_specialq` builds all `(p, roots)` progressions for the coefficient's
`[p_min, p_max]` range into `p_soa_array_t` (SoA by root-count, degree-4..7
have 2/4/8, 1/5/25, etc. — see `p_soa_array_init`), giving `t->num_entries`
total (p, root) pairs. It then pulls special-q from `sieve_fb_next` in
batches.

### 2.2 Batch sizing — memory-bound only

```c
max_batch_specialq32 = d->max_sort_entries32 / t->num_entries;
max_batch_specialq64 = d->max_sort_entries64 / t->num_entries;
```

`max_sort_entries{32,64}` come from `gpu_sieve_data_init`, sized from
`0.3 * gpu_mem` (or a `gpu_mem_mb=` override) divided across `num_threads`,
capped at 50M/35M elements. **This cap is purely a GPU-memory/throughput
figure — it has no dependence on `key_bits`.** That mismatch is the root
cause in §7.

Per-batch, `max_batch_size` is additionally clipped to `1 << unused_bits`
(spare high bits in the 32-bit `p_max` word, used to disambiguate which
`aprog` value a hit came from).

### 2.3 `key_bits` derivation

```c
key_bits = ceil(log((double)p_max * p_max * ((num_aprog_vals + 1) / 2)) / M_LN2);
if (num_aprog_vals > 1) key_bits++;
```

`key_bits` depends only on `p_max` (fixed for the whole `sieve_specialq`
call) and `num_aprog_vals` (normally 1; boosted only on a short tail batch
where `batch_size < max_batch_size / 3`, to keep the GPU fed by having each
`(p,q)` pair emit multiple offsets). **`key_bits` never looks at `batch_size`
in the dominant `num_aprog_vals == 1` case** — it's essentially a constant
for the whole coefficient, set by the target's bit-size (smaller factorization
targets ⇒ smaller `p_max` ⇒ smaller `key_bits`).

`root_bytes` (32 vs 64-bit keys on the wire) is decided once from `key_bits
> 32`.

### 2.4 Handoff (`handle_special_q_batch`)

Launches the `sieve_kernel_trans_*` kernel to materialize roots into
`t->gpu_root_array`, then either:
- `d->use_collision_engine`: fills `collision_data_t` (`num_elements =
  num_specialq * t->num_entries * num_aprog_vals`, `key_bits`, `root_bytes`,
  `bucket_hash`, `debug`/`collect_stats` flags) and calls
  `collision_engine_run`, or
- else: the plain CUB-sort engine path (`sort_engine_run` +
  `sieve_kernel_final_*`) — no bucketing, no `CANDIDATE_CAP`, not subject to
  §7's failure mode.

Collision stats (`bucket_max`, `candidate_count`, `dedup_count`,
`value_match_count`, grow/hash-cap counts, `filter_iters_hist[102]`,
`elapsed_ms`) accumulate into `device_thread_data_t` and get logged at the
end of `sieve_specialq` when `collision_stats` is on (histogram broken into
zero/converged/cap-hit buckets, plus a bucket-size-vs-iteration-count
breakdown).

Engine selection/config comes from `nfs_args` tokens (`read_collision_engine_args`):
`collhash=`, `collstats=`, `colldebug=`, `colllib=` (override the `.so`/`.dll`
path). The GPU library requires compute capability 7.0+ (Volta); missing lib
load prints that requirement explicitly.

---

## 3. Device side: the five-stage kernel pipeline (`collision_engine.cu`)

`collision_engine_run(engine, data)` — one call per special-q batch:

```
scatter_roots_kernel        bucket every root by a hash of the key
        │                   (retry loop: grow max_per_bucket on overflow)
        ▼
filter_per_bucket_kernel    per-bucket, iteratively narrow toward true
        │                   duplicates via widening hash slices (§4)
        │                   → writes surviving items to candidate_keys,
        │                     atomicAdd's their count into candidate_cnt
        ▼
cub::DeviceRadixSort        sort all candidates globally
        ▼
dedup_kernel                keep only keys with a neighbor duplicate
        │                   (is_first && has_neighbor) → dedup_keys
        ▼
count_secondary_kernel      build a second-level hash table (D/S) over
scatter_secondary_kernel    the deduped keys → X (bucketed by hash_value)
        ▼
count_matched_values_*      scan the ORIGINAL n roots again, look each up
scatter_matched_values_*    in X via find_candidate_slot, collect all
  (or the _arena variant)   (p, q) values sharing a colliding key
        ▼
emit_found_kernel           for each colliding key-slot, pairwise-check
  (or _arena variant)       gcd(p1,p2)==1 and q1==q2, store_hit_collision
```

### 3.1 Bucket scatter (`scatter_roots_kernel`)

Hashes each root (`compute_bucket`, from `collision_bucket.h`) into one of
`NUM_BUCKETS = 1 << LOG2_NUM_BUCKETS` buckets. Two scatter modes:
`SCATTER_DIRECT_ATOMIC` (one atomicAdd per thread) or the default warp-
aggregated path (`__match_any_sync`/`__ballot`-style leader election, one
atomicAdd per distinct bucket-value per warp — cheaper under bucket
contention). Slots beyond `max_per_bucket` set `overflow_flag` and drop the
item (recovered by the retry loop below, not lost silently).

### 3.2 Bucket-capacity retry loop (`collision_engine_run`)

```c
for (;;) {
    scatter_roots_kernel<<<...>>>(...);
    cub::DeviceReduce::Max(..., d_bucket_count, d_max_bucket, ...);
    if (!overflow && max_bucket <= max_per_bucket) break;
    grown_cap = observed + observed/4 + 64;
    engine->ensure_capacity(n, key_bits, grown_cap);   // reallocs + rescatters
}
```

Self-correcting: if the actual max bucket occupancy exceeds the current
allocation, it grows 25%+64 past the observed max and rescatters from
scratch. `data->bucket_grow_count` tracks how often this fires.

### 3.3 Per-bucket iterative filter (`filter_per_bucket_kernel`) — the core

One CUDA block per bucket, dynamic shared memory holds three bitmap tables
`T`, `T2`, `T3` (each `max_tsize_words` words), sized by
`compute_capped_ilog2(cnt, key_bits, max_tsize_words)` — i.e. table size
scales with observed bucket occupancy but is capped by both remaining key
bits (`key_bits - LOG2_NUM_BUCKETS`) and available shared memory
(`cudaDevAttrMaxSharedMemoryPerBlockOptin`, ≤ ~⅓ of opt-in max split three
ways).

**Pre-pass** (before the iteration loop): hash every item in the bucket on
its lowest `ilog2` bits above `LOG2_NUM_BUCKETS` (`my_shift2 =
LOG2_NUM_BUCKETS` initially) into `T2` (seen-once) / `T` (seen-twice, i.e.
"this hash value has a same-bucket collision").

**Iteration loop** (`MAX_FILTER_ITERS = 20`): each round
1. re-hashes only the *surviving* items from last round, checking
   `U[hv>>5] & bit` (was this item's *previous* hash value flagged as
   colliding) — this is the actual filter step;
2. items that survive get re-hashed on the *next* `ilog2`-bit slice
   (`my_shift2 = my_shift + 6`) into a *new* pair of tables, ping-ponged
   between `{T,T3}` on even/odd iterations so the "read" table from last
   round and the "write" table for next round never alias;
3. warp-ballot compaction (`__ballot_sync` + `__popc` prefix) packs
   survivors into `arr_out` without per-thread atomics on the common path.

**Stop conditions** (`stop_zero`, `stop_cap` at `it == MAX_FILTER_ITERS-1`,
or `stop_conv` when `s_nsize[it-3] == s_nsize[it]`): whichever fires first,
the current survivor set is emitted as "candidates" via a single
`atomicAdd(candidate_cnt, cnt)` **before** checking capacity — so
`candidate_cnt` accumulates the true post-filter survivor count from *every*
bucket regardless of whether the global cap is later found to be exceeded
(see §7).

**Shift wraparound**: `my_shift2 + ilog2 > key_bits` resets `my_shift2` back
to `LOG2_NUM_BUCKETS`, i.e. re-hashing on bits already consumed by the
pre-pass. Once this triggers, further iterations add no new discriminating
information — `stop_conv` will fire on a stabilized-but-not-actually-small
survivor count. This is expected/harmless when the bucket is small enough
that a handful of rounds already isolated true duplicates before wraparound;
it's the mechanism behind §7 when it isn't.

Also tracked (`iters_hist_out`, 102 slots): stop-iteration histogram (total /
zero / cap-hit, 21 bins each) and a 13-bin log2 bucket-size histogram split
by how many iterations it took (fast ≤3 / medium ==4 / slow ≥5) — this is
what `collision_stats` prints as `filter iters converged:` / `filter
bucket-size ...:` in the host log.

### 3.4 Global dedup

Candidates from all buckets are radix-sorted (`cub::DeviceRadixSort`,
straight 64-bit keys), then `dedup_kernel` keeps exactly the keys with a
non-unique neighbor in sorted order (`is_first && has_neighbor`) — this is
where survivor noise from an unresolved filter (§3.3, §7) turns into a much
larger `dedup_cnt`/downstream cost than a healthy run would show.

### 3.5 Secondary hash + value gathering

`count_secondary_kernel`/`scatter_secondary_kernel` build a classic
counting-sort index (`D` = per-slot counts → exclusive-scanned to offsets,
`S` = occupied-slot bitmap, `X` = the deduped keys themselves, ordered by
slot) over the (typically much smaller) `dedup_cnt` keys — sized by
`c_ilog2 = host_ilog2(dedup_cnt+1)`, capped at `MAX_C_ILOG2 = 20`.

Then the **original** `n`-element root/value stream is scanned again
(`count_matched_values_kernel` / `scatter_matched_values_kernel`, or the
fused `count_and_store_matched_values_kernel` "arena" variant) — for each
root, `find_candidate_slot` binary-searches `X[D[hv]..D[hv+1])` for an exact
key match, and on hit records the associated `(p,q)` value.

**Arena vs. general path**: if `dedup_cnt <= VALUE_MATCH_CAP /
MATCH_ARENA_WIDTH`, a fixed-width-8 arena per slot is tried first
(`count_and_store_matched_values_kernel`, single pass, no exclusive-scan
needed) — cheaper, but overflows (an actual duplicate hash value with >8
matches) fall back to the general two-pass path (count → exclusive-sum →
scatter). `data->match_arena_fallback_count` tracks how often that
fallback fires.

### 3.6 Emit

`emit_found_kernel` (general path, uses `value_offsets`) / `emit_found_arena_kernel`
(arena path, uses fixed-width slots) do the final pairwise check per colliding
slot: split each matched value into `(q, p)` via `pshift`, require `q1==q2`
and `gcd(p1,p2)==1`, and on success call `store_hit_collision` — an
atomicAdd-allocated slot in the shared `found_t` ring buffer (index 0 is a
running counter; capacity `found_array_size`, silently drops past-capacity
hits rather than corrupting memory — host side, `check_found_array` in
`stage1_sieve_gpu.c` reads this back and reports saturation via
`found_saturated_batches`).

---

## 4. Compile-time constants (`collision_engine.cu`)

| constant             | value                    | role                                       |
|----------------------|--------------------------|---------------------------------------------|
| `MAX_FILTER_ITERS`   | 20                       | per-bucket filter round cap                |
| `BLOCK_THREADS`      | 128                      | threads/block for scatter + filter kernels |
| `CANDIDATE_CAP`      | `1u << 22` = 4,194,304   | global cap on post-filter survivors — **fixed, independent of `n`/`key_bits`** |
| `VALUE_MATCH_CAP`    | = `CANDIDATE_CAP`        | cap on matched `(p,q)` values              |
| `MATCH_ARENA_WIDTH`  | 8                        | per-slot fixed capacity in the arena fast path |
| `MAX_C_ILOG2`        | 20                       | cap on secondary hash-table `ilog2`        |
| `LOG2_NUM_BUCKETS`   | 14 (`collision_bucket.h`)| `NUM_BUCKETS = 16384` — fixed regardless of `key_bits` |

`ensure_capacity(n, key_bits, min_per_bucket)` reallocs everything sized off
`max_n`/`max_key_bits`/`max_per_bucket` only when a new call needs more than
what's already allocated (monotonic growth, never shrinks within a run).
Initial `max_per_bucket` estimate: `mean + 6*sqrt(mean) + 32` (mean =
`n / NUM_BUCKETS`), i.e. ~6-sigma of a Poisson/normal approximation before
the retry-and-grow loop (§3.2) kicks in for skewed distributions.

---

## 5. `collision_data_t` field cheat-sheet (host ↔ device contract)

Set by `handle_special_q_batch`, consumed by `collision_engine_run`:

| field           | meaning                                                     |
|-----------------|--------------------------------------------------------------|
| `keys_in`       | `t->gpu_root_array` — the generated roots (32 or 64-bit)     |
| `data_in`       | `t->gpu_p_array` — the `(p_index)`/value stream, same length |
| `q_batch`       | device special-q array, indexed by the `q` extracted from a matched value |
| `found_array`   | device `found_t` ring buffer                                 |
| `num_elements`  | `num_specialq * t->num_entries * num_aprog_vals` = `n`       |
| `key_bits`      | bit-width of the key domain (§2.3) — **the crux of §7**      |
| `root_bytes`    | 4 or 8, selects 32/64-bit key path throughout                |
| `shift`         | `32 - unused_bits`, used by `emit_found_*` to split matched values into `(q,p)` |
| `bucket_hash`   | hash-mode selector passed to `compute_bucket`                |
| `debug` / `collect_stats` | gate `printf` diagnostics / stats accumulation     |

---

## 6. Where things live

| file                                                        | role |
|--------------------------------------------------------------|------|
| `gnfs/poly/stage1/stage1_sieve_gpu.c`                        | host driver: batching, `key_bits`, thread/context lifecycle, engine dynamic-load, stats logging |
| `gnfs/poly/stage1/stage1_core_gpu/collision_engine.cu`        | device kernels + `collision_engine` class; built into `cub/collision_engine.{so,dll}`, loaded at runtime (requires CC 7.0+) |
| `collision_bucket.h`                                          | `LOG2_NUM_BUCKETS`, `NUM_BUCKETS`, `BUCKET_MASK`, `compute_bucket()` — shared with any future bucket-scatter kernel |
| `stage1_core_gpu/stage1_core.h`                               | `found_t`, `specialq_t`, `collision_data_t` struct defs (the host/device contract) |
| `FUSED_TRANS_SCATTER_PLAN.md`                                 | retired fused trans+scatter spike (2026-05-25, reverted — register pressure regression) — historical, not live code |

---

## 7. Known failure mode: candidate overflow on small inputs

**Symptom**: `collision_engine: candidate overflow <N> > 4194304` and
`exit(-1)` in `collision_engine_run`, with `<N>` a large multiple of
`CANDIDATE_CAP` (observed: 49,683,408 ≈ 11.8×, at `n=49,999,884`,
`key_bits=23`, `LOG2_NUM_BUCKETS=14`). Triggers mainly on small inputs
(< ~480 bits).

**Root cause, corrected**: my first-pass diagnosis (a clean `n / 2^key_bits`
occupancy-ratio story, with a `key_bits`-scaled batch cap as the fix) was an
oversimplification — the actual relationship between special-q batch size,
`key_bits`, `shift` (`32 - unused_bits`, the unused top-bit count of
`p_max`), and the resulting `n` fed to the collision engine is more
entangled than that (batch size interacts with `key_bits` and `shift`
together, and `n` can grow disproportionately fast as batch size grows —
not a clean linear/ratio relationship). Not fully re-derived; treat the
clean formula above as **wrong**, not just approximate.

**Fix actually applied**: rather than compute a `key_bits`-derived cap, just
hard-cap the special-q batch size outright when the collision engine is in
use, in `sieve_specialq` where `max_batch_size` is set (§2.2):

```c
if (d->use_collision_engine) {
	// hard cap for gpu_gerbicz, which has some built-in
	// caps on candidate counts that too-large of batch
	// can exceed.  triggers mainly on small inputs (< 480 bits).
	max_batch_size = MIN(max_batch_size, 16384);
}
```

This resolved the overflow in practice. `16384` was arrived at empirically,
not derived from `CANDIDATE_CAP`/`key_bits`/`shift` analytically — the exact
relationship between batch size and worst-case `n` (and thus how much
headroom this constant actually has, or whether it needs to vary with
`key_bits`/`shift` for other target sizes) is still an open question worth
revisiting if overflows resurface at a different size regime.

The GPU-underutilization risk flagged in earlier drafts of this doc is a
non-issue in practice: this cap only bites on small inputs (< ~480 bits),
where poly hits are plentiful and cheap regardless — losing some GPU
efficiency there doesn't cost wall-clock time the way it would on a large,
hit-starved target. That's part of why the empirical hard cap is an
acceptable stopping point rather than something that needs the clean
derivation.

---

### The one-line mental model

Roots get scattered into 16384 fixed hash buckets; each bucket iteratively
narrows toward its true duplicates using ever-more-significant key bits
(ping-ponged bitmap tables, stopping on empty/cap/no-progress); survivors
funnel through a global sort+dedup+second-hash to reattach their `(p,q)`
values and get pairwise gcd/q-checked into `found_t`. Everything scales with
`n` (memory-bound); nothing scales with `key_bits` — which is fine when
`key_bits` is large (sparse collisions) and breaks (§7) when it isn't.
