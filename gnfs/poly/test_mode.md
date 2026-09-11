# Stage-1 test mode — forced-window hit tap for engine cross-check

A flag-guarded diagnostic that makes every engine sweep an **identical** forced
`(a_d, p-window, q-window)` and dump its raw valid hits `(a_d, p, m)` to a
per-engine file, **without running stage 2**. Off by default; when off, the code
path is byte-identical to production. No engine internals, factory, sort, or
CUB/Gerbicz code is touched — the tap sits in the shared `handle_collision`, and
the window is forced through the *existing* `specialq` range parameters.

Run once per engine with the same tokens, then `sort -u | comm` the dumps.

## Tokens (via `nfs_stage1_args`)

```
test_ad=<decimal a_d>     # presence of this enables test mode
test_pmin=<uint32>  test_pmax=<uint32>
test_qmin=<uint64>  test_qmax=<uint64>
```
Pick a window inside the common envelope (p < 2^27, q < 2^32 — the ground-truth
tier from the fixtures) so no engine's `cell_fits` clamps it; then all four get
the identical window.

---

## 1. `stage1.h` — test config on `stage1_sieve_data_t`

```c
	/* stage-1 cross-check test mode (0 = off) */
	int    test_mode;
	mpz_t  test_ad;
	uint32 test_pmin, test_pmax;
	uint64 test_qmin, test_qmax;
	FILE  *test_dump;
```

## 2. `stage1.c` — parse + arm, in `stage1_sieve_data_init`

After the engine is selected (`v` known) and `num_threads` is computed, before
`d->num_threads = num_threads;` and the pool creation:
```c
	d->test_mode = 0;
	mpz_init(d->test_ad);
	d->test_dump = NULL;
	if (obj->nfs_args != NULL) {
		const char *tmp;
		if ((tmp = strstr(obj->nfs_args, "test_ad=")) != NULL) {
			char b[160]; uint32 n = 0;
			tmp += 8;
			while (n < sizeof(b)-1 && *tmp && *tmp != ' ' && *tmp != ',')
				b[n++] = *tmp++;
			b[n] = 0;
			mpz_set_str(d->test_ad, b, 10);
			d->test_mode = 1;
		}
		if ((tmp = strstr(obj->nfs_args, "test_pmin=")) != NULL) d->test_pmin = strtoul (tmp+10, NULL, 10);
		if ((tmp = strstr(obj->nfs_args, "test_pmax=")) != NULL) d->test_pmax = strtoul (tmp+10, NULL, 10);
		if ((tmp = strstr(obj->nfs_args, "test_qmin=")) != NULL) d->test_qmin = strtoull(tmp+10, NULL, 10);
		if ((tmp = strstr(obj->nfs_args, "test_qmax=")) != NULL) d->test_qmax = strtoull(tmp+10, NULL, 10);
	}
	if (d->test_mode) {
		char fn[256];
		num_threads = 1;                         /* single writer -> no dump lock */
		sprintf(fn, "test_%s.hits", v->name);
		d->test_dump = fopen(fn, "w");
		logprintf(obj, "TEST MODE %s: a_d=%Zd  p[%u,%u]  q[%" PRIu64
			",%" PRIu64 "]  -> %s\n", v->name, d->test_ad,
			d->test_pmin, d->test_pmax, d->test_qmin, d->test_qmax, fn);
	}
```

## 3. `stage1.c` — single-a_d bypass in `search_coeffs`

Right after `poly_coeff_t *c = poly_coeff_init();` and `deadline_per_coeff = …;`,
**before** `init_ad_sieve(...)` (so it skips the whole a_d generator + pre-count):
```c
	if (d->test_mode) {
		mpz_set(c->high_coeff, d->test_ad);
		stage1_bounds_update(poly, c);    /* sets m0/coeff_max/trans_* for handle_collision */
		search_coeff_async(d, c, deadline_per_coeff);
		poly_coeff_free(c);
		return;                           /* the one task drains in data_free */
	}
```

## 4. `stage1.c` — force the window in `search_coeff_core`

**(a)** right after `special_q_min = 1;`:
```c
	if (d->test_mode) {
		p_min = d->test_pmin;
		p_max = d->test_pmax;
		special_q_min = d->test_qmin;
		special_q_max = d->test_qmax;
	}
```
**(b)** keep it one piece — change the `num_pieces` guard:
```c
	if (!d->test_mode && special_q_max - special_q_min > 500000)   /* +!d->test_mode */
		num_pieces = MIN(200, ...);
```
Now `sieve_fb_init` builds each engine's factory for the test q-range, `num_pieces`
stays 1, and the dispatch calls `specialq` with exactly `[test_qmin,test_qmax] ×
[test_pmin,test_pmax]`. `cell_fits` still runs (leave it) — with an in-envelope
window it won't clamp, so every engine gets the same range.

## 5. `stage1.c` — the tap in `handle_collision`

Right after the final `mpz_tdiv_q(c->m, c->m, c->tmp1);` (m fully reconstructed),
**before** the `{ /* submit … */ }` block:
```c
	if (task->d->test_mode) {
		if (task->d->test_dump)
			gmp_fprintf(task->d->test_dump, "%Zd %Zd %Zd\n",
					c->high_coeff, c->p, c->m);
		return;                           /* skip stage 2 entirely */
	}
```
This dumps only *valid* hits (past every Kleinjung check) as `(a_d, p, m)` — the
same triple the old `.m` files used, so the `comm` workflow is unchanged.

## 6. `stage1.c` — close the dump in `stage1_sieve_data_free`

After the pools are drained/freed (so all writes are flushed), before the struct
is torn down:
```c
	if (d->test_dump)
		fclose(d->test_dump);
	mpz_clear(d->test_ad);
```

---

## Running the cross-check

Four runs, identical window, one per engine:
```
… nfs_stage1_args="stage1_engine=cpu_gerbicz   test_ad=<A> test_pmin=<P0> test_pmax=<P1> test_qmin=<Q0> test_qmax=<Q1>"
… nfs_stage1_args="stage1_engine=cpu_hashtable  test_ad=<A> …same window…"
… nfs_stage1_args="stage1_engine=gpu_cubsort    test_ad=<A> …same window…"
… nfs_stage1_args="stage1_engine=gpu_gerbicz    test_ad=<A> …same window…"
```
Each drops `test_<engine>.hits`. Then:
```
for f in cpu_gerbicz cpu_hashtable gpu_cubsort gpu_gerbicz; do sort -u test_$f.hits > $f.s; done
comm -3 cpu_gerbicz.s cpu_hashtable.s | wc -l   # same CPU factory -> expect ~0
comm -3 gpu_cubsort.s  gpu_gerbicz.s   | wc -l   # same GPU factory -> expect ~0
comm -3 cpu_gerbicz.s  gpu_cubsort.s   | head    # CPU vs GPU -> the divergence, quantified
comm -23 cpu_gerbicz.s gpu_cubsort.s | wc -l     # in CPU only
comm -13 cpu_gerbicz.s gpu_cubsort.s | wc -l     # in GPU only
```

## Reading it

- **within a family ~identical** (`cpu_gerbicz`≈`cpu_hashtable`,
  `gpu_cubsort`≈`gpu_gerbicz`): the engines are correct — same factory, same
  window, same hits. Any non-trivial diff *inside* a family is a real engine bug.
- **across families differs**: expected — the CPU engines use the driver's
  factory (`sieve_fb_init` at 100–5000 / 2–200000), the GPU engines use their own
  in `gpu_thread_data_init`. This mode **quantifies** that divergence directly:
  same window in, how different are the hits out. That number is what tells you
  whether the §7 factory-normalization is worth doing (it is if the families
  disagree) and gives you a before/after target once you do it.

Note: test mode forces one thread and never touches stage 2, so it's fast and
deterministic — but it's a *diagnostic*, not a benchmark. Timing means nothing
here; only the hit sets do.
