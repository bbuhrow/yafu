# Wiring `poly_stats` into the poly-select pipeline

Adds one shared stats instance, reached from the driver (hits) and the
stage-2 callbacks (sizeopt / rootopt / best-E), plus queue depth. Build
`poly_stats.c` with the tree (see §6). All hooks are guarded `if (…stats)`,
so an unset pointer is a no-op.

Anchors are surrounding code, not line numbers. `#1`
(`stage1_sieve_cpu_hashtable.c`) is already patched separately.

---

## 1. `gnfs/poly/poly_skew.h`

**1a — include the module** (after the existing `poly.h` include):
```c
#include "poly.h"
#include "poly_stats.h"          /* ADD */
```

**1b — `poly_stage1_t`: carry the stats pointer to the driver**
```c
	stage1_callback_t callback;
	void *callback_data;
	poly_stage_stats_t *stats;    /* ADD */
} poly_stage1_t;
```

**1c — `rootopt_callback_data_t`**
```c
typedef struct {
	FILE* all_poly_file;
	poly_config_t* config;
	poly_stage_stats_t *stats;    /* ADD */
} rootopt_callback_data_t;
```

**1d — `sizeopt_callback_data_t`**
```c
typedef struct {
	poly_rootopt_t* rootopt;
	rootopt_callback_data_t* rootopt_callback;
	poly_stage_stats_t *stats;    /* ADD */
} sizeopt_callback_data_t;
```

---

## 2. `gnfs/poly/poly_skew.c`

**2a — `rootopt_callback`: count a saved poly + best-E.** Insert between the
`fflush(data->all_poly_file);` and `save_poly(config, &poly);`:
```c
	fflush(data->all_poly_file);

	if (data->stats) {                                   /* ADD */
		char adbuf[64];                              /* ADD */
		gmp_snprintf(adbuf, sizeof(adbuf), "%Zd",    /* ADD */
				coeff1[degree]);             /* ADD */
		poly_stats_add_rootopt(data->stats,          /* ADD */
				combined_score, adbuf);      /* ADD */
	}                                                    /* ADD */

	save_poly(config, &poly);
```

**2b — `sizeopt_callback`: count a size-opt survivor.**
```c
	sizeopt_callback_data_t *callback = (sizeopt_callback_data_t *)extra;

	if (callback->stats)                                 /* ADD */
		poly_stats_add_sizeopt(callback->stats);     /* ADD */

	poly_rootopt_run(callback->rootopt, alg_coeffs,
			rat_coeffs, sizeopt_norm, projective_alpha);
```

**2c — orchestrator (the `poly_skew`/`find_poly_skew` function).**

Declarations, beside the other locals (`… rootopt_callback_data_t
rootopt_callback_data;`):
```c
	rootopt_callback_data_t rootopt_callback_data;
	poly_stage_stats_t poly_stats;    /* ADD */
	int poly_verbose = 0;             /* ADD */
```

Parse verbosity, in the `if (obj->nfs_args != NULL)` block (after the
`max_coeff=` parse):
```c
		tmp = strstr(obj->nfs_args, "poly_verbose=");  /* ADD */
		if (tmp != NULL)                               /* ADD */
			poly_verbose = atoi(tmp + 13);         /* ADD */
```

Wire + report — replace the lone `poly_stage1_run(obj, &stage1_data);`:
```c
		poly_stats_init(&poly_stats, poly_verbose, 0);          /* ADD */
		stage1_data.stats = &poly_stats;                        /* ADD */
		if (obj->flags & MSIEVE_FLAG_NFS_POLYSIZE)              /* ADD */
			sizeopt_callback_data.stats = &poly_stats;      /* ADD */
		if (obj->flags & MSIEVE_FLAG_NFS_POLYROOT)              /* ADD */
			rootopt_callback_data.stats = &poly_stats;      /* ADD */

		poly_stage1_run(obj, &stage1_data);                     /* existing */

		poly_stats_report(&poly_stats, 1);                      /* ADD */
		poly_stats_free(&poly_stats);                           /* ADD */
```
(The counting all happens inside `poly_stage1_run`, so the whole lifecycle
sits in the `MSIEVE_FLAG_NFS_POLY1` block.)

---

## 3. `gnfs/poly/stage1/stage1.h`  (polysize)

**`stage1_sieve_data_t`** — beside the `engine` field added earlier:
```c
	const struct stage1_engine_vtable *engine;
	poly_stage_stats_t *stats;        /* ADD */
} stage1_sieve_data_t;
```
(`poly_stage_stats_t` is visible via `poly_skew.h` → `poly_stats.h`, which
`stage1.h` already includes transitively.)

---

## 4. `gnfs/poly/stage1/stage1.c`  (polysize driver)

**4a — `stage1_hit_data_t`: carry stats to the stage-2 worker**
```c
typedef struct {
	stage1_callback_t callback;
	void *callback_data;
	poly_stage_stats_t *stats;        /* ADD */

	mpz_t ad;
	mpz_t p;
	mpz_t m;
} stage1_hit_data_t;
```

**4b — `stage1_hit_run`: mark task complete (drives queue depth)**
```c
	hit_data->callback(hit_data->ad,
			   hit_data->p,
			   hit_data->m,
			   hit_data->callback_data);

	if (hit_data->stats)                                 /* ADD */
		poly_stats_stage2_done(hit_data->stats);     /* ADD */
}
```

**4c — `handle_collision`: count the hit + hand stats to the task.** In the
submit block:
```c
		hit_data->callback = d->poly->callback;
		hit_data->callback_data = d->poly->callback_data;
		hit_data->stats = d->stats;                  /* ADD */
		mpz_init_set(hit_data->ad, c->high_coeff);
		mpz_init_set(hit_data->p, c->p);
		mpz_init_set(hit_data->m, c->m);

		if (d->stats)                                /* ADD */
			poly_stats_add_hit(d->stats);        /* ADD */

		threadpool_add_task(d->stage2_threadpool,
					&task_control, 1);
```

**4d — `poly_stage1_run`: attach the instance to the driver data**
```c
	stage1_sieve_data_init(&sieve_data, obj, &poly);
	sieve_data.stats = data->stats;                      /* ADD */

	search_coeffs(&sieve_data, data->deadline);
```

**4e — `search_coeffs`: pre-count a_d (for the fraction) + mark current a_d.**

Add a counter beside the other locals:
```c
	double cumulative_time = 0;
	uint64 ad_index = 0;              /* ADD */
```

Pre-count, inserted **after** the existing begin-alignment block and
**before** `while (1) {`. It walks a scratch sieve, then restores the search
state, so the real loop is byte-identical to before — worst case the
denominator is off by one at the range boundary, never the search itself:
```c
	mpz_mul_ui(poly->gmp_high_coeff_begin, poly->tmp1,
			ad_sieve.high_coeff_multiplier);

	/* ADD: count the a_d in range for the progress fraction */
	if (d->stats) {
		sieve_t count_sieve;
		poly_coeff_t *cc = poly_coeff_init();
		mpz_t save_begin;
		uint64 ad_count = 0;

		mpz_init_set(save_begin, poly->gmp_high_coeff_begin);
		init_ad_sieve(&count_sieve, poly);
		while (find_next_ad(&count_sieve, poly, cc->high_coeff) == 0)
			ad_count++;
		free_ad_sieve(&count_sieve);
		poly_coeff_free(cc);
		mpz_set(poly->gmp_high_coeff_begin, save_begin);
		mpz_clear(save_begin);
		d->stats->ad_total = ad_count;
	}
	/* END ADD */

	while (1) {
```

Mark the current a_d, right after `find_next_ad` succeeds:
```c
		if (find_next_ad(&ad_sieve, poly, c->high_coeff))
			break;

		if (d->stats) {                                      /* ADD */
			char adbuf[64];                              /* ADD */
			gmp_snprintf(adbuf, sizeof(adbuf), "%Zd",    /* ADD */
					c->high_coeff);              /* ADD */
			poly_stats_set_ad(d->stats, adbuf, ++ad_index); /* ADD */
		}                                                    /* ADD */

		stage1_bounds_update(poly, c);
```

---

## 5. `factor/nfs/nfs_poly.c`  (YAFU) — map `-v`/`-v -v` to `poly_verbose`

In the `#else` (non-CUDA) branch, after the `nfs_args` composition and the
`nfs_stage1_args` append you already added, append the verbosity token from
YAFU's own flag:
```c
	sprintf(nfs_args + strlen(nfs_args), " poly_verbose=%d",
			t->fobj->VFLAG);
```
Now `-v` → roll-up, `-v -v` → scrolling roll-up, no flag → silent — the same
`VFLAG` YAFU already uses everywhere.

---

## 6. Build (`Makefile`)

Add the object beside `poly_skew`, in the same object list / with the same
suffix the poly objects use (e.g. the `.no` variant list `stage1_engine`
went into):
```
$(BUILDDIR)/gnfs/poly/poly_stats$(OBJ_EXT)
```
`poly_stats.c` compiles with `-std=gnu*` like the rest (it uses no inline asm,
but keep the dialect consistent).

---

## Sanity after applying

- `-v` prints one overwriting line; `-v -v` scrolls it; no `-v` is silent.
- The `q` (queue-depth) field should sit low at first and, if it climbs and
  stays high, that's your stage-2-bound signal — the thing the upcoming
  parallel stage 2 is meant to drain.
- At end of run `q` returns to 0 (every submitted hit's task completed).
- `hits ≥ sizeopt ≥ rootopt` always (a strict funnel); if not, a hook is
  double-firing.
