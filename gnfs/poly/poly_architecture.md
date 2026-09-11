# NFS Polynomial Selection — Architecture Reference

A concise map of the current (post-refactor) poly-select pipeline: how work
flows, how threads are laid out, how the callback chain nests, and which knobs
control the CPU behaviour. Written for the YAFU + msieve-polysize integration.

---

## 1. What it does, in one breath

Poly select finds the best degree-d algebraic polynomial for an NFS
factorization. It runs in two phases per leading coefficient `a_d`:

- **Stage 1 (collision search)** — find `(p, m)` pairs that make a promising
  raw polynomial. Cheap-ish, high volume, embarrassingly parallel over `a_d`.
- **Stage 2 (optimization)** — for each raw hit, **size-optimize** it, then
  **root-optimize** it (a root sieve over rotations), scoring each result by
  Murphy-E and saving the good ones. Expensive, and where the wall time goes.

Stage 1 **produces** hits; stage 2 **consumes** them. A bounded blocking queue
sits between, so the whole thing is a producer/consumer pipeline.

---

## 2. Threading topology (one instance owns everything)

```
  YAFU  do_msieve_polyselect
    │   forces fobj->THREADS = 1  (one YAFU worker over the whole a_d range)
    │   passes the REAL thread count to msieve as num_threads
    ▼
  msieve  find_poly_core  →  poly_stage1_run  →  stage1_sieve_data_init
    │
    ├─ stage1_threadpool : num_threads workers   ← PRODUCERS (one a_d each)
    │        │  each runs the selected collision engine over its a_d
    │        │  every collision → handle_collision → submit
    │        ▼
    │   ┌──────────────────────────────┐
    │   │ bounded blocking queue (1000)│   ← back-pressure: producers BLOCK
    │   └──────────────────────────────┘      when full  (this is the balancer)
    │        │
    ├─ stage2_threadpool : S workers    ← CONSUMERS (size-opt + root-opt)
    │        S = stage2_threads token (default 1)
    ▼
  status line + one shared .p file
```

Key point: it's **one** msieve instance now (the override collapsed the old
"THREADS YAFU workers each spawning a THREADS-wide pool" nesting). All threads
— `num_threads` producers + `S` consumers — live in this single instance, so
the stats, best-E, and `.p` file are genuinely job-wide.

The blocking queue is itself a load balancer: when stage 2 is behind, producers
block on the full queue and yield their cores; when it's caught up, consumers
block on the empty queue. Runnable-thread count tracks the busy side.

---

## 3. End-to-end flow

### Producer side (stage 1)
```
poly_stage1_run
 └ search_coeffs                         (stage1.c)  — the a_d loop
    └ find_next_ad                        → next smooth leading coeff a_d
    └ stage1_bounds_update                → norm_max ⇒ coeff_max, p_size_max,
    │                                        m0, sieve_size   (the stage-1 gate)
    └ search_coeff_async → (task) search_coeff_core
         └ stage1_engine_cell_fits        → skip/clamp to the engine envelope
         └ v->specialq(...)               → the selected engine's collision search
              └ (per collision) handle_collision
                   ├ Kleinjung modular consistency check, build (a_d, p, m)
                   ├ poly_stats_add_hit                     (hits++)
                   └ threadpool_add_task(stage2, blocking)  → into the queue
```

### Consumer side (stage 2)
```
stage1_hit_run(data, threadid)           (stage1.c)  — a stage-2 worker
 └ extra = num_stage2_workers>0 ? &stage2_workers[threadid].sizeopt_data
                                : hit_data->callback_data   (single-worker path)
 └ stage1_callback(ad, p, m, extra)      → poly_sizeopt_run(sizeopt_data, …)
      └ (per size-opt candidate ≤ max_sizeopt_norm)  sizeopt_callback
           ├ poly_stats_add_sizeopt                  (sizeopt++)
           └ poly_rootopt_run(rootopt_data, …)       → the root sieve
                └ (per rotation with combined_E > min_e)  rootopt_callback
                     ├ [file_lock] write "# norm … e …" block to .p; save_poly
                     └ poly_stats_add_rootopt         (rootopt++, best-E)
 └ poly_stats_stage2_done                (stage2_done++  → drives queue depth)
```

**The funnel is not monotone.** `hits ≥ cand` (most collisions die at
size-opt), but root-opt **fans out**: one size-opt candidate spawns a whole
root sieve that saves *many* polynomials, so `saved ≥ cand` is normal.
`saved` = the line count in the `.p` file.

---

## 4. The callback nest (why it looks confusing)

Each phase is written as "do work, then hand each result to a callback you were
given at init." Three phases ⇒ three registered callbacks, chained. The
indirection is that a phase never names the next phase directly — it just calls
`data->callback`.

```
 registration (in find_poly_core / build_stage2_workers)
 ─────────────────────────────────────────────────────────
 poly_stage1_init (&stage1_data,  stage1_callback,  &sizeopt_data)
 poly_sizeopt_init(&sizeopt_data, sizeopt_callback, &sizeopt_callback_data)
 poly_rootopt_init(&rootopt_data, rootopt_callback, &rootopt_callback_data)

 invocation chain (per hit)
 ─────────────────────────────────────────────────────────
 stage1_callback(extra = sizeopt_data)
     → poly_sizeopt_run(sizeopt_data)
         → sizeopt_data.callback == sizeopt_callback(extra = sizeopt_callback_data)
             → poly_rootopt_run(sizeopt_callback_data.rootopt == rootopt_data)
                 → rootopt_data.callback == rootopt_callback(extra = rootopt_callback_data)
                     → write .p
```

What each `extra` (callback_data) actually carries — this is the part worth
keeping straight:

| callback          | its `extra` is            | which holds                                              |
|-------------------|---------------------------|---------------------------------------------------------|
| `stage1_callback` | `poly_sizeopt_t*`         | the size-opt working state (private per stage-2 worker) |
| `sizeopt_callback`| `sizeopt_callback_data_t*`| `→ rootopt_data`, `→ rootopt_callback_data`, `stats`    |
| `rootopt_callback`| `rootopt_callback_data_t*`| `all_poly_file`, `config`, `file_lock`, `stats`         |

Parallel stage 2 works by giving **each consumer its own
`stage2_worker_t`** = `{ sizeopt_data, rootopt_data, sizeopt_callback_data,
rootopt_callback_data }`. The optimization state is private; only
`all_poly_file`, the poly heap (`save_poly`), best-E, and stats are shared —
serialized by `file_lock` (writes) and the stats mutex (counters). The router
in `stage1_hit_run` picks `stage2_workers[threadid]` by the pool-local thread
id, so no two consumers touch the same buffers.

---

## 5. Stage-1 engines (runtime-selectable)

Four collision engines sit behind one vtable contract, chosen by the
`stage1_engine=` token; `HAVE_CUDA` gates *compilation* of the GPU pair, the
registry gates *selection* at run time.

| engine          | id | device | notes                                    |
|-----------------|----|--------|------------------------------------------|
| `cpu_gerbicz`   | #2 | CPU    | default; sort-based; the reference       |
| `cpu_hashtable` | #1 | CPU    | trunk hashtable, wrapped (`HAVE_CPU_HASHTABLE`) |
| `gpu_cubsort`   | #3 | GPU    | CUB radix sort (`HAVE_CUDA`)             |
| `gpu_gerbicz`   | #4 | GPU    | kyleaskine collision engine (`HAVE_CUDA`)|

Each engine supplies `*_thread_data_init/free` and
`stage1_specialq_<engine>(task, threadid, q_min, q_max, p_min, p_max)`, all
feeding the one shared `handle_collision`. An **envelope** (`max_p`,
`max_special_q`) per engine lets the driver skip an over-cap `a_d` or clamp an
over-cap `q` window before dispatch. (Empirically, on CPU the engine choice
barely moves end-to-end yield — stage 2 dominates — which is the argument for
GPU stage 1: move the producer off the CPU so all cores serve stage 2.)

---

## 6. Instrumentation

`poly_stage_stats_t` (one shared, mutex-guarded instance) tracks the funnel and
prints a single job-wide roll-up line:

```
a_d 372 (17/68) | hits 1.24M  cand 47.1k  saved 312 | q 5 | best E 1.281e-9 @ a_d 804 | 84s
```

- `hits / cand (sizeopt candidates passing stage2 norm check) / saved (rootopt passing min Murphy-E check)` — the three funnel counters (§3).
- `q = hits − stage2_done` — the stage-2 backlog (queued + in-flight). Pegged
  near 1000 ⇒ stage 2 can't keep up; raise `stage2_threads`.
- verbosity from `poly_verbose` (mapped from YAFU `VFLAG`): 0 silent · 1 in-place
  roll-up · 2 scrolling roll-up.

---

## 7. CPU / tuning knobs

All the `*=` tokens ride the `nfs_stage1_args` passthrough (set in `yafu.ini`
or on the CLI); they land in `obj->nfs_args` and are parsed where noted.

| knob                | set via                    | controls                                            |
|---------------------|----------------------------|-----------------------------------------------------|
| `-t N`              | YAFU CLI                   | stage-1 producer count (`num_threads`)              |
| `stage2_threads=S`  | `nfs_stage1_args`          | stage-2 consumer count (`build_stage2_workers`)     |
| `stage1_engine=…`   | `nfs_stage1_args`          | collision engine (registry)                         |
| `poly_verbose=N`    | auto from `VFLAG`          | roll-up verbosity                                   |
| `stage1_norm=…`     | `nfs_args`                 | stage-1 gate (search window, `stage1_bounds_update`)|
| `stage2_norm=…`     | `nfs_args`                 | size-opt gate (`pol_norm·e^alpha ≤ max_sizeopt_norm`)|

Balancing rule of thumb: total useful threads ≈ cores, split producers vs
consumers. Because root-opt fans out (heavy consumer side), the sweet spot
usually has `S ≥ num_threads`. The blocking queue self-limits, so
over-provisioning thread *count* doesn't runaway — a blocked worker yields its
core. Tune by walking `stage2_threads` up until `q` stops pegging without CPU
overshooting your core count.

---

## 8. The four run modes

`find_poly_core` branches on the msieve flags (YAFU `-np*`):

| mode        | flags                    | behaviour                                             |
|-------------|--------------------------|------------------------------------------------------|
| `-np`       | POLY1+POLYSIZE+POLYROOT  | full pipeline; **parallel stage 2 engages** at `S>1` |
| `-np1`      | POLY1 only               | stage 1 → dump hits to `.m`; one stage-2 writer       |
| `-nps`      | POLYSIZE only            | sequential loop: read `.m`, size-opt each             |
| `-npr`      | POLYROOT only            | sequential loop: read size-opt file, root-opt each    |

Only `-np` with `stage2_threads>1` takes the concurrent per-worker-bundle path.
The standalone modes are sequential file readers using a single
sizeopt/rootopt instance — which is why `rootopt_callback` must no-op its
`file_lock`/`stats` when they're `NULL` (the `memset` of the callback-data
structs guarantees that).

---

## 9. Where things live

| file                                   | role                                              |
|----------------------------------------|---------------------------------------------------|
| `factor/nfs/nfs_poly.c`                | YAFU driver; the `THREADS=1` override             |
| `gnfs/poly/poly_skew.c`                | `find_poly_core`; the three callbacks; `build_stage2_workers` |
| `gnfs/poly/stage1/stage1.c`            | search driver: `search_coeffs`, `stage1_bounds_update`, `handle_collision`, `stage1_hit_run`, the pools |
| `gnfs/poly/stage1/stage1_engine.{c,h}` | engine registry + envelopes                       |
| `gnfs/poly/stage1/stage1_sieve_cpu.c`  | #2 `cpu_gerbicz` (reference)                       |
| `gnfs/poly/stage1/stage1_sieve_cpu_hashtable.c` | #1 `cpu_hashtable` (wrapped trunk)       |
| `gnfs/poly/poly_stats.{c,h}`           | shared stats + roll-up reporter                   |
| `gnfs/poly/stage2/optimize.c`          | size-opt + Murphy scoring (`rootopt_callback` fires here) |
| `gnfs/poly/stage2/stage2.c`            | root-opt entry (`sizeopt_callback` fires here)    |
| `gnfs/poly/stage2/root_sieve*.c`       | the root sieve (the fan-out)                      |

---

## 10. Key structs at a glance

- `poly_stage1_t` — stage-1 interface: N, degree, `norm_max`, `deadline`,
  `callback`+`callback_data`, `stage2_workers`+`num_stage2_workers`, `stats`.
- `stage1_sieve_data_t` (`d`) — the running search: `obj`, `poly`, `num_threads`,
  `threads[]`, the two pools, `engine`, `stage2_workers`, `stats`.
- `task_data_t` — one `a_d`'s search task (`obj`, `c`, `d`, deadline).
- `stage1_hit_data_t` — one submitted hit (`callback`, `callback_data`, `d`,
  `stats`, `ad/p/m`).
- `stage2_worker_t` — one consumer's private bundle (§4).
- `poly_stage_stats_t` — the shared funnel + reporter (§6).

---

### The one-line mental model
One instance; `num_threads` producers walk `a_d` and fire collisions through
`handle_collision` into a blocking queue; `S` consumers each own a private
size-opt/root-opt bundle and drain the queue, fanning out into the shared `.p`
under a lock; `q` is the tension between the two, and `stage2_threads` is the
dial. GPU work will change who the producers are, not this shape.
