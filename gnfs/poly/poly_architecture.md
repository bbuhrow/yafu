# NFS Polynomial Selection — Architecture Reference

A concise map of the current poly-select pipeline: how work flows, how
threads are laid out, how the callback chain nests, and which knobs control
CPU vs. GPU behaviour. Written for the YAFU + msieve-polysize integration.

**Update note:** stage 1 → stage 2 dispatch is no longer uniform across
engines. GPU keeps the original queue/pool design; CPU now runs stage 2
inline on the hit-finding thread. See §2–§4 for the split.

---

## 1. What it does, in one breath

Poly select finds the best degree-d algebraic polynomial for an NFS
factorization. It runs in two phases per leading coefficient `a_d`:

- **Stage 1 (collision search)** — find `(p, m)` pairs that make a promising
  raw polynomial. Cheap-ish, high volume, embarrassingly parallel over `a_d`.
- **Stage 2 (optimization)** — for each raw hit, **size-optimize** it, then
  **root-optimize** it (a root sieve over rotations), scoring each result by
  Murphy-E and saving the good ones. Expensive, and where the wall time goes.

Stage 1 **produces** hits; stage 2 **consumes** them. How that handoff
happens now differs by engine:

- **GPU** — a bounded blocking queue sits between producer and consumers,
  same as before: a genuine producer/consumer pipeline.
- **CPU** — there is no queue. Each producer thread runs stage 2 itself,
  synchronously, the moment it finds a hit, then resumes its own collision
  search. "Producer" and "consumer" are the same thread.

---

## 2. Threading topology

```
  YAFU  do_msieve_polyselect
    │   forces fobj->THREADS = 1  (one YAFU worker over the whole a_d range)
    │   passes the REAL thread count to msieve as num_threads
    ▼
  msieve  find_poly_core  →  poly_stage1_run  →  stage1_sieve_data_init
    │
    ├─ stage1_threadpool : num_threads workers   ← one a_d each
    │        │  each runs the selected collision engine, in q-window chunks
    │        │  every collision → handle_collision
    │        ▼
    │   ┌────────────────────── branch on engine ──────────────────────┐
    │   │                                                               │
    │   │  CPU engine                         GPU engine                │
    │   │  ──────────                         ──────────                │
    │   │  run stage 2 inline,                submit to stage2_threadpool│
    │   │  same thread, right now             (blocking queue, size 1000)│
    │   │  using this producer's own          │                         │
    │   │  stage2_workers[threadid]            ▼                         │
    │   │  bundle, then resume                stage2_threadpool : S      │
    │   │  collision search                   workers  ← CONSUMERS      │
    │   │                                      (size-opt + root-opt)     │
    │   └───────────────────────────────────────────────────────────────┘
    │        │
    ▼
  status line + one shared .p file
```

Still **one** msieve instance owning everything — stats, best-E, and the
`.p` file are job-wide regardless of engine. What changed is only the
stage-1→stage-2 handoff mechanism, not the shared state around it.

**GPU** keeps the original load-balancing property: producers block on a
full queue and yield their cores; consumers block on an empty queue.

**CPU** has no such balancing to do — there's no idle thread to move work
to, since a thread that isn't searching is by definition busy running
stage 2 for the hit it just found. The `q` backlog metric (§6) is
consequently always ~0 for CPU-only runs; it only means something for GPU.

---

## 3. End-to-end flow

### Producer side (stage 1) — same shape for both engines
```
poly_stage1_run
 └ search_coeffs                         (stage1.c)  — the a_d loop
    └ find_next_ad                        → next smooth leading coeff a_d
    └ stage1_bounds_update                → norm_max ⇒ coeff_max, p_size_max,
    │                                        m0, sieve_size   (the stage-1 gate)
    └ search_coeff_async → (task) search_coeff_core
         └ stage1_engine_cell_fits        → skip/clamp to the engine envelope
         └ [CPU only] chunk [q_min,q_max) into STAGE1_QCHUNK_DIVISOR pieces;
              after each chunk, check wall-clock deadline
              (task->coeff_deadline) and abort the a_d early if exceeded
         └ v->specialq(...)               → the selected engine's collision search,
                                              once per chunk (CPU) or once total (GPU)
              └ (per collision) handle_collision(task, threadid, ...)
                   ├ Kleinjung modular consistency check, build (a_d, p, m)
                   ├ poly_stats_add_hit                     (hits++)
                   └ dispatch — branches by engine (see below)
```

**Dispatch inside `handle_collision`, by engine:**

```c
if (d->engine->envelope.is_gpu) {
    /* build hit_data, threadpool_add_task(d->stage2_threadpool, ..., 1) */
}
else {
    /* CPU: call d->poly->callback(...) directly, right here, using
       &d->stage2_workers[threadid].sizeopt_data; then
       poly_stats_stage2_done(d->stats) immediately, no queue involved */
}
```

### Consumer side (stage 2)

**GPU** — unchanged from before:
```
stage1_hit_run(data, threadid)           (stage1.c)  — a stage-2 worker
 └ extra = &stage2_workers[threadid].sizeopt_data   (threadid = consumer id)
 └ stage1_callback(ad, p, m, extra)      → poly_sizeopt_run(sizeopt_data, …)
      └ (per size-opt candidate ≤ max_sizeopt_norm)  sizeopt_callback
           ├ poly_stats_add_sizeopt                  (sizeopt++)
           └ poly_rootopt_run(rootopt_data, …)       → the root sieve
                └ (per rotation with combined_E > min_e)  rootopt_callback
                     ├ [file_lock] write "# norm … e …" block to .p; save_poly
                     └ poly_stats_add_rootopt         (rootopt++, best-E)
 └ poly_stats_stage2_done                (stage2_done++)
```

**CPU** — the same callback chain, but reached directly from
`handle_collision` (no `stage1_hit_run`, no queue, no separate consumer
thread — `threadid` here is the *producer's own* id):
```
handle_collision (CPU branch)
 └ extra = &stage2_workers[threadid].sizeopt_data   (threadid = producer id)
 └ d->poly->callback(ad, p, m, extra)     → poly_sizeopt_run(sizeopt_data, …)
      └ ... identical chain to the GPU path from here down ...
 └ poly_stats_stage2_done                (called inline, right after)
```

**The funnel is not monotone**, same as before. `hits ≥ cand` (most
collisions die at size-opt), but root-opt **fans out**: one size-opt
candidate spawns a whole root sieve that saves *many* polynomials, so
`saved ≥ cand` is normal. `saved` = the line count in the `.p` file. This
is unaffected by the dispatch change — same callback chain either way.

---

## 4. The callback nest (why it looks confusing)

Unchanged in shape from before — three phases, three registered callbacks,
chained; a phase never names the next phase directly, it calls
`data->callback`:

```
 registration (in find_poly_core / build_stage2_workers)
 ─────────────────────────────────────────────────────────
 poly_stage1_init (&stage1_data,  stage1_callback,  &sizeopt_data)
 poly_sizeopt_init(&sizeopt_data, sizeopt_callback, &sizeopt_callback_data)
 poly_rootopt_init(&rootopt_data, rootopt_callback, &rootopt_callback_data)

 invocation chain (per hit) — identical for CPU and GPU from here down
 ─────────────────────────────────────────────────────────
 stage1_callback(extra = sizeopt_data)
     → poly_sizeopt_run(sizeopt_data)
         → sizeopt_data.callback == sizeopt_callback(extra = sizeopt_callback_data)
             → poly_rootopt_run(sizeopt_callback_data.rootopt == rootopt_data)
                 → rootopt_data.callback == rootopt_callback(extra = rootopt_callback_data)
                     → write .p
```

What each `extra` (callback_data) actually carries — unchanged:

| callback          | its `extra` is            | which holds                                              |
|-------------------|---------------------------|-----------------------------------------------------------|
| `stage1_callback` | `poly_sizeopt_t*`         | the size-opt working state (private per `stage2_worker_t`) |
| `sizeopt_callback`| `sizeopt_callback_data_t*`| `→ rootopt_data`, `→ rootopt_callback_data`, `stats`      |
| `rootopt_callback`| `rootopt_callback_data_t*`| `all_poly_file`, `config`, `file_lock`, `stats`           |

**What changed: who `stage2_workers[]` is indexed by, and how big it is.**
`stage2_worker_t` (private `sizeopt_data`/`rootopt_data`/callback-data;
shared `all_poly_file`, `file_lock`, `stats`) is now used by **both**
GPU consumer threads and CPU producer threads, from the same array — just
never both at once in a single run (a run uses one engine). The array is
sized `MAX(stage2_threads, num_threads)` (previously just `stage2_threads`)
so there's always one bundle per CPU producer thread even when
`stage2_threads=` was never set. Index meaning depends on which engine is
active:
- **GPU:** index = consumer thread id (from `stage2_threadpool`), as before.
- **CPU:** index = producer thread id (from `stage1_threadpool`) — the
  same thread that found the hit uses its own bundle directly, no routing
  needed since there's only ever one candidate index per thread.

---

## 5. Stage-1 engines (runtime-selectable)

Four collision engines sit behind one vtable contract, chosen by the
`stage1_engine=` token; `HAVE_CUDA` gates *compilation* of the GPU pair, the
registry gates *selection* at run time.

| engine          | id | device | notes                                    | stage-2 dispatch |
|-----------------|----|--------|-------------------------------------------|------------------|
| `cpu_gerbicz`   | #2 | CPU    | default; sort-based; the reference       | inline, same thread |
| `cpu_hashtable` | #1 | CPU    | trunk hashtable, wrapped (`HAVE_CPU_HASHTABLE`) | inline, same thread |
| `gpu_cubsort`   | #3 | GPU    | CUB radix sort (`HAVE_CUDA`)             | queued to `stage2_threadpool` |
| `gpu_gerbicz`   | #4 | GPU    | kyleaskine collision engine (`HAVE_CUDA`)| queued to `stage2_threadpool` |

Each engine supplies `*_thread_data_init/free` and
`stage1_specialq_<engine>(task, threadid, q_min, q_max, p_min, p_max)`, all
feeding the one shared `handle_collision` (now also taking `threadid`). An
**envelope** (`max_p`, `max_special_q`, `is_gpu`) per engine lets the driver
skip an over-cap `a_d`, clamp an over-cap `q` window, and — new — decide
`handle_collision`'s dispatch branch. CPU engines additionally have their
`v->specialq()` call chunked by `search_coeff_core` (§3, §7); GPU keeps one
whole-window call, since chunking is a deadline-enforcement mechanism for
inline CPU work and GPU's own runtime already isn't blocked on stage 2.

---

## 6. Instrumentation

`poly_stage_stats_t` (one shared, mutex-guarded instance) tracks the funnel and
prints a single job-wide roll-up line:

```
a_d 372 (17/68) | hits 1.24M  cand 47.1k  saved 312 | q 5 | best E 1.281e-9 @ a_d 804 | 84s
```

- `hits / cand (sizeopt candidates passing stage2 norm check) / saved (rootopt passing min Murphy-E check)` — the three funnel counters (§3), meaningful identically for both engines.
- `q = hits − stage2_done` — the stage-2 backlog (queued + in-flight).
  **GPU:** pegged near 1000 ⇒ stage 2 can't keep up; raise `stage2_threads`.
  **CPU:** `poly_stats_stage2_done` fires immediately after the inline
  callback returns, so `q` stays ~0 by construction — there's no backlog to
  read a signal from on the CPU path; it's not a tuning lever there.
- verbosity from `poly_verbose` (mapped from YAFU `VFLAG`): 0 silent · 1 in-place
  roll-up · 2 scrolling roll-up.

---

## 7. CPU / tuning knobs

All the `*=` tokens ride the `nfs_stage1_args` passthrough (set in `yafu.ini`
or on the CLI); they land in `obj->nfs_args` and are parsed where noted.

| knob                | set via                    | controls                                            | engine relevance |
|---------------------|----------------------------|------------------------------------------------------|-------------------|
| `-t N`              | YAFU CLI                   | producer count (`num_threads`)                        | both |
| `stage2_threads=S`  | `nfs_stage1_args`          | GPU consumer-pool size (`build_stage2_workers`)       | GPU only — CPU always gets one bundle per producer thread regardless of this setting |
| `stage1_engine=…`   | `nfs_stage1_args`          | collision engine (registry)                           | both |
| `coeff_deadline=…`  | `nfs_stage1_args`          | per-`a_d` wall-clock time budget, checked at chunk boundaries | CPU only — GPU doesn't chunk |
| `poly_verbose=N`    | auto from `VFLAG`          | roll-up verbosity                                     | both |
| `stage1_norm=…`     | `nfs_args`                 | stage-1 gate (search window, `stage1_bounds_update`)  | both |
| `stage2_norm=…`     | `nfs_args`                 | size-opt gate (`pol_norm·e^alpha ≤ max_sizeopt_norm`) | both |

**Balancing rule of thumb — now engine-dependent:**
- **GPU:** unchanged — total useful threads ≈ cores, split producers vs
  consumers; because root-opt fans out, the sweet spot usually has
  `S ≥ num_threads`. The blocking queue self-limits, so over-provisioning
  thread *count* doesn't run away. Tune by walking `stage2_threads` up
  until `q` stops pegging without CPU overshooting your core count.
- **CPU:** no such tuning exists — every producer thread is always also
  its own stage-2 consumer, so `num_threads` alone determines both. The
  relevant tuning is `coeff_deadline` (bound one `a_d`'s search+optimize
  time) and the internal `STAGE1_QCHUNK_DIVISOR`/`STAGE1_QCHUNK_MAX`
  chunk-size constants (not yet exposed as `nfs_args` tokens), which trade
  off deadline-check granularity against per-chunk call overhead.

---

## 8. The four run modes

`find_poly_core` branches on the msieve flags (YAFU `-np*`):

| mode        | flags                    | behaviour                                             |
|-------------|--------------------------|--------------------------------------------------------|
| `-np`       | POLY1+POLYSIZE+POLYROOT  | full pipeline; CPU dispatches stage 2 inline per hit; GPU still needs `S>1` for the concurrent per-worker-bundle path |
| `-np1`      | POLY1 only               | stage 1 → dump hits to `.m`; no stage 2 at all, so the CPU/GPU dispatch split doesn't apply here |
| `-nps`      | POLYSIZE only            | sequential loop: read `.m`, size-opt each — unaffected, single sizeopt instance regardless of engine |
| `-npr`      | POLYROOT only            | sequential loop: read size-opt file, root-opt each — unaffected |

Only `-np` is affected by the CPU/GPU dispatch split, since it's the only
mode where stage 1 and stage 2 run concurrently in the same process. The
standalone modes are unchanged sequential file readers using a single
sizeopt/rootopt instance — `rootopt_callback` still no-ops its
`file_lock`/`stats` when they're `NULL` for those paths.

---

## 9. Where things live

| file                                   | role                                              |
|----------------------------------------|---------------------------------------------------|
| `factor/nfs/nfs_poly.c`                | YAFU driver; the `THREADS=1` override             |
| `gnfs/poly/poly_skew.c`                | `find_poly_core`; the three callbacks; `build_stage2_workers` (now sized `MAX(stage2_threads, num_threads)`) |
| `gnfs/poly/stage1/stage1.c`            | search driver: `search_coeffs`, `stage1_bounds_update`, `search_coeff_core` (chunking, CPU only), `handle_collision` (engine dispatch branch), `stage1_hit_run` (GPU path only now), the pools |
| `gnfs/poly/stage1/stage1_engine.{c,h}` | engine registry + envelopes (`is_gpu` now also drives dispatch, not just device residency) |
| `gnfs/poly/stage1/stage1_sieve_cpu.c`  | #2 `cpu_gerbicz`; `finish_search` is the actual `handle_collision` call site; internal per-1000-root deadline check (needs wall-clock, not `get_cpu_time()`) |
| `gnfs/poly/stage1/stage1_sieve_cpu_hashtable.c` | #1 `cpu_hashtable` (wrapped trunk)       |
| `common/util.c`                        | `get_cpu_time()` (process-wide, existing) and `get_wall_time()` (new — `CLOCK_MONOTONIC`, used for all deadline checks) |
| `gnfs/poly/poly_stats.{c,h}`           | shared stats + roll-up reporter                   |
| `gnfs/poly/stage2/optimize.c`          | size-opt + Murphy scoring (`rootopt_callback` fires here) |
| `gnfs/poly/stage2/stage2.c`            | root-opt entry (`sizeopt_callback` fires here)    |
| `gnfs/poly/stage2/root_sieve*.c`       | the root sieve (the fan-out)                      |

---

## 10. Key structs at a glance

- `poly_stage1_t` — stage-1 interface: N, degree, `norm_max`, `deadline`,
  `callback`+`callback_data`, `stage2_workers`+`num_stage2_workers`, `stats`.
- `stage1_sieve_data_t` (`d`) — the running search: `obj`, `poly`,
  `num_threads`, `threads[]`, both pools (`stage2_threadpool` idle/unused
  on CPU-engine runs), `engine`, `stage2_workers` (now dual-purpose — see
  §4), `stats`.
- `task_data_t` — one `a_d`'s search task (`obj`, `c`, `d`, `coeff_deadline`
  — now actually enforced via wall-clock at chunk boundaries, not dead).
- `stage1_hit_data_t` — one submitted hit (`callback`, `callback_data`, `d`,
  `stats`, `ad/p/m`) — **GPU path only** now; CPU no longer constructs these.
- `stage2_worker_t` — one bundle's private sizeopt/rootopt state — now
  shared in concept between a GPU consumer and a CPU producer, indexed by
  whichever role is active (§4).
- `poly_stage_stats_t` — the shared funnel + reporter (§6); `q` only
  carries balancing information for GPU.

---

### The one-line mental model
One instance; `num_threads` producers walk `a_d` in wall-clock-bounded
chunks. On a **collision**, GPU engines queue the hit for a separate
consumer pool (`S` workers, tuned via `stage2_threads`, balanced by a
blocking queue); CPU engines instead run the entire stage-2 chain
immediately, on the same thread, using that thread's own private bundle,
then resume searching. Either way the callback chain, the shared `.p` file
under `file_lock`, and the stats are identical from `stage1_callback`
onward — only how a hit gets there, and on which thread, differs by engine.
