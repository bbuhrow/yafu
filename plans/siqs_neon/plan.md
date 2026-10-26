# Plan: AArch64 NEON kernels for SIQS

*Started 2026-10-25 on branch `L2sieve` (after `f26dee4`). Target machine:
Neoverse-N1 (AWS Graviton2, 2 cores, 64 KB L1d, 1 MB L2/core, 32 MB L3,
NEON + dotprod, no SVE).*

## 0. Status

| phase | item | status |
|---|---|---|
| P0 | local cross build + qemu test loop | done: `tools/neon_cross_build.sh`, `tools/neon_ab.sh` |
| P0 | baseline profile on the N1 | done (0.2): med sieve 46%, root updates 15% |
| P0 | NEON on/off switch for A/B checks | done: `-noneon`, `fact_obj_t.HAS_NEON` |
| P0 | `neon_util.h` (movemask, mulhi_u16) | done |
| P0 | ARM builds without `FORCE_GENERIC` | done (see 0.1) |
| P1 | relation scanner | done (`check_relations_siqs_neon`); qemu A/B matches; N1 timing pending |
| P2 | nextRoots (poly root updates) | small primes done (`nextRoots_32k_neon_small`, exact under qemu); mid buckets dropped: NEON now defaults to the sparse mid sieve |
| P3 | med sieve | done: `split` is the default (C73: 39.6 s vs 41.1 s `ranges`, 48.2 s `-noneon`); med-sieve-only gain pending (`YAFU_NEON_MED=portable`) |
| P4 | tdiv_med + resieve_med | not started |
| P5 | tdiv_LP (bucket tdiv) | not started |
| P6 | firstRoots pre-pass, lp_sieveblock | only if the profile says so |
| P7 | Lanczos `v_t` ops | last |
| P8 | tuning: params table, 64k blocks | after the kernels |

What already works on ARM (as of `f26dee4`): the whole build, SIQS with
correct factors, the 4-lane direct-sieve large-prime kernels
(`direct_sieve_neon.c` via `direct_sieve_v4.h`), and CPU/cache detection.
Everything else runs the portable C paths. The banner shows this as
"Generic C enabled 32k sieve core" and "portable root updates".

Reference point: C70 `rsa(232)`, `-noopt`, 1 thread: 27.5 s sieving
(`FORCE_GENERIC=1` build).

### 0.1 Build notes

- The first N1 builds used `FORCE_GENERIC=1` because `MOD_CMP_8X` (tdiv
  med) and `spPRP2` exist only as x86 asm. Now `USE_8X_MOD_ASM` and
  `SSE2_RESIEVING` are x86-only (qs_impl.h), and tdiv.c uses the GMP
  PRP check where `spPRP2` is missing, so a plain `make` works on ARM.
  Resieving (`USE_RESIEVING`, portable path) stays on. On a C61 the
  plain build and the `FORCE_GENERIC` build give identical statistics.
- `common/vec_bitonic_sort.c` is x86-64 only (only the KNL root updates
  use it).
- **Always A/B with `-noopt`.** The SPV tuner (SIQS.c, `tf_small_cutoff`)
  picks its value from measured rels/sec, so without `-noopt` the same
  seed can give different polys/reports/relations (seen often under
  qemu).
- `-noneon` reverts every NEON kernel (scanner, direct-sieve LP kernels,
  and later ones) to the portable code, so `tools/neon_ab.sh` compares
  the whole NEON set against portable in one binary.

### 0.2 N1 baseline profile (2026-10-25)

Plain build (no `FORCE_GENERIC`), with the NEON scanner (1.78%),
`siqs(rsa(240))` `-threads 1 -seed 1 -noopt`, `perf record` on cycles,
top of `--no-children --sort symbol`:

| share | symbol | NEON status |
|---|---|---|
| 45.98% | `med_sieveblock_32k` | portable C -> **P3, do first** |
| 12.17% | `direct_resieve_update_creg_neon` | NEON already (half is `direct_compress_flush`) |
| 11.03% | `nextRoots_32k_lp_buckets` | portable -> P2 |
| 4.10% | `nextRoots_32k_generic_small` | portable -> P2 |
| 3.91% | `firstRoots_prep_generic` | portable scalar pre-pass -> P6 |
| 3.82% | `direct_resieve_hits` | scalar |

The med sieve is a much bigger share than on x86 (33% with AVX-512BW),
so the phase order below is now P3, P2, then the rest by profile.

## 1. Approach

- **Port from the intrinsic versions, not the inline asm.** Most x86 SIQS
  kernels exist twice: GCC inline asm (`GCC_ASM64X`) and an MSVC
  intrinsic build of the same algorithm. The intrinsic versions (and the
  SSE4.1 files) are 128-bit and mostly 16-bit lanes, which map almost
  1:1 onto NEON (`vaddq_u16`, `vqsubq_u16`, `vceqq_u16`, `vminq_u16`, ...).
- **Selection.** NEON is mandatory on AArch64, so there is no runtime
  check. Add `*_neon` kernels behind
  `#if defined(__aarch64__) && defined(__ARM_NEON)` and assign them in
  the dispatch block of `siqs_static_init` (SIQS.c, around the
  `med_sieve_ptr`/`nextRoots_ptr` assignments) and next to the
  `scan_ptr` choice.
- **Keep the portable path selectable** (e.g. `-siqsLPisa generic` or a
  new `-noneon`) so every kernel can be A/B'd on the same seed.
- **Shared helpers** go in one header (`factor/qs/neon_util.h`):
  movemask equivalents (`vshrn_n_u16` to a 64-bit nibble mask, or
  `vaddvq` with bit weights as `direct_sieve_neon.c` does), the unsigned
  16x16 high multiply (`vmull_u16` + `vuzp2q_u16`, since NEON has no
  `pmulhuw`), and compress via `vqtbl1q_u8`.
- **No x86-only flags leak in**: kernels compile with plain
  `-march=armv8-a`; dotprod only behind `__ARM_FEATURE_DOTPROD`.

Missing instructions to plan around:

| x86 | NEON replacement |
|---|---|
| `pmovmskb` / `movemask_epi8` | `vshrn_n_u16(vreinterpretq_u16_u8(m), 4)` gives a 64-bit mask, 4 bits per byte |
| `pmulhuw` | `vmull_u16` (lo/hi halves) then `vuzp2q_u16` |
| `vpcompressd` | 16-entry `vqtbl1q_u8` LUT (already used) |
| `tzcnt` / `blsr` | `__builtin_ctzll` (`rbit`+`clz`) / `x & (x-1)` |
| scatter / gather | scalar |

## 2. Inventory: x86 SIMD kernels vs. the current ARM path

| stage | x86 kernels | ARM now | NEON port source |
|---|---|---|---|
| relation scan | `check_relations_siqs_16_{sse2,avx2,avx512}` (tdiv_scan.c) | C fallback of `_sse2` | SSE2/AVX2 scan: OR blocks, test high bits, then find bytes |
| poly root update, med primes | `nextRoots_32k_sse41`, `_avx2(_intrin)` | `nextRoots_32k` scalar | `poly_macros_common_sse4.1.h` + MSVC paths in update_poly_roots_32k.c |
| poly root update, bucket primes | AVX2/AVX-512 bucket fill | scalar | AVX2 intrin version, 4 lanes of 32 bits |
| first roots (per A) | `firstRoots_vec_avx512/avx2`, scalar pre-pass | scalar pre-pass | the scalar pre-pass; maybe `float64x2_t` modmul later |
| med sieve | `med_sieveblock_32k_{sse41,avx2,avx512bw}(_split)` | `med_sieveblock_32k` C | `med_sieve_32k_sse4.1.c`, `sieve_macros_32k_sse4.1.h` (MSVC branch) |
| tdiv med primes | `tdiv_medprimes_32k(_avx2)`, SSE2 blocks under `D_HAS_SSE2` | C | SSE2 16-bit Barrett check: needs the `pmulhuw` replacement |
| resieve med primes | `resieve_medprimes_32k_{avx2,avx512bw}`, SSE2 | C | same as tdiv med |
| tdiv large (bucket) | `tdiv_LP_{sse2,avx2,avx512}` (tdiv_large.c) | C branch of `_sse2` | AVX2 compare of bucket fb/loc fields, 4 lanes |
| bucket sieve dump | `lp_sieveblock_avx512{f,bw}` | `lp_sieveblock` C | scatter-bound; probably stays scalar |
| small-prime tdiv (SPV) | AVX-512BW under `TRY_SMALLP_OPT` | C | only if profiled hot |
| direct large primes | creg/cmem kernels | **NEON (done)** | n/a |
| Lanczos | msieve `v_t` ops, `VBITS=128` | C on 2x64 | `uint64x2_t` `veorq`/`vandq`; small share at SIQS sizes |

## 3. Testing

1. **Local correctness loop (no cloud round trip).** In WSL (Ubuntu
   22.04) with `gcc-aarch64-linux-gnu` and `qemu-user-static`:
   `tools/neon_cross_build.sh` syncs the checkout to `~/neon/yafu` and
   cross-builds a static aarch64 yafu (GMP 6.2.1 cross-built to
   `~/neon/gmp-arm64`); `tools/neon_ab.sh ~/neon/yafu 'siqs(rsa(200))'`
   runs NEON vs `-noneon` and compares statistics. A C61 takes ~10 s
   under qemu. qemu timings mean nothing; the same `neon_ab.sh` runs
   natively on the N1.
2. **Bit-exact kernel checks.** Each new kernel gets a verify build
   (as `FIRSTROOTS_VERIFY` does): run NEON and portable on the same input
   and compare outputs (roots, sieve bytes, reports, relations).
   `SIQS_NEON_EMU` + `arm_neon_emu.h` stay available for x86-side checks
   of small kernels.
3. **Performance on the N1 only.** `perf stat -e cycles:u` and
   cycles/poly with a fixed `-seed`, 1 thread, on c70/c75/c80, since wall
   time is noisy. If the VM has no PMU, use `-e cpu-clock` and more
   repetitions.

## 4. Phases

Order after the N1 profile (0.2): P3 (med sieve, 46%) first, then P2
(root updates, 15%), then P4-P6 by their measured shares.

- **P0 infrastructure:** cross toolchain + qemu loop; the NEON on/off
  switch; `neon_util.h`; the baseline `perf record` on the N1.
- **P1 scanner:** small, self-contained, and it settles the movemask idiom
  that the later kernels reuse. Done: `check_relations_siqs_neon`
  (tdiv_scan.c) replaces the 4/8/16 variants for every size: max-reduce
  of 64 bytes, exact movemask only on a hit. Same report rules as the
  SSE2 scanner (skip loc 0 and >= 65534); with TLP/QLP it compares
  against `tf_closnuf` like the AVX512BW scanner. The scanner was a
  small share on x86, so expect a small gain.
- **P2 nextRoots:** med primes from the SSE4.1 macros, then bucket primes
  from the AVX2 intrinsic code. The x86 survey suggests this plus the med
  sieve is most of the non-direct time. Status:
  - `nextRoots_32k_neon_small` (update_poly_roots_32k_neon.c) replaces
    `nextRoots_32k_generic_small` (4.7%) on the `lp_portable_roots` path:
    `[sieve_small_fb_start, med_B)` in 8 lanes of 16 bits; the sign < 0
    update is done as r - (p - u) + borrow, since r + u can exceed 16
    bits. Exact-logp builds still match `-noneon` on the C61 and C70.
  - The 12% `nextRoots_32k_lp_buckets` was the mid-range bucket fill,
    used only because the `wide` test in `siqs_resolve_strategy` left
    NEON out (buckets by default). On the N1 (C73) `-siqsMid sparse` was
    42.0 s vs 43.7 s for buckets, so NEON is now in `wide`: the sparse
    mid sieve (NEON kernels, vector sparse root update) is the default
    and no NEON bucket fill is needed. `-noneon` (scalar ISA) still gets
    buckets, as before.
  - A NEON bucket fill anyway, in case it beats sparse:
    `nextRoots_32k_lp_buckets_neon` (update_poly_roots_32k_neon.c),
    called from `nextRoots_32k_lp_neon` via `V4_LP_BUCKETS` (SSE4.1 and
    the NEON emulation keep the portable fill). Roots 4 primes at a time;
    the fill goes prime by prime after each `CHECK_NEW_SLICE`, so slices,
    bucket contents and element order are identical to the portable fill.
    Fill modes (`YAFU_NEON_BUCKETS`): `fixed` (default) runs
    K = interval / p_min + 1 branch-free passes per root stream, a miss
    writing bucket 0's next free slot without counting it; `loop` is the
    portable while loop. Exact builds match `-noneon` with
    `-siqsMid bucket` in both modes.
- **P3 med sieve:** from the SSE4.1 file; then try the fixed-pass `_split`
  structure that won on AVX2 and AVX-512. Done, in
  `factor/qs/med_sieve_32k_neon.c`:
  - `med_sieveblock_32k_neon` ("ranges"): the SSE4.1 structure. Scalar
    below 2^13, then 8 lanes of 16-bit roots with 3/2/1/0 unconditional
    steps + 1 conditional for 8192-10922, -2^14, -2^15, >2^15. The
    `fb_*_B` bounds are the first multiple-of-8 index *past* each
    threshold, so each vector loop checks its limit on the group's last
    prime (`full_fb->list->prime`, unaffected by a-prime zeroing); the
    first version trusted the bounds and gave 2 extra reports on a C61.
  - `med_sieveblock_32k_neon_split`: the AVX2 `_split` structure. 8-lane
    dump/selective/last-root1 loops from start_prime to `fb_32k_div3`
    (so primes < 8192 are vectorized too), then branch-free fixed passes
    (K = 3/2/1 by the group's smallest prime) with the sink strip at
    `sieve[-256, 0)`.
  - logp: one value per range (ranges: 13/14/15) or per 8-group (split:
    the group's first prime), like x86. `-DMED_NEON_EXACT_LOGP` uses each
    prime's own logp; with it both variants match `-noneon` exactly on
    the C61 and C70 under qemu.
  - Selection: `MED_SIEVE_NEON_FN` (default `med_sieveblock_32k_neon_split`),
    overridden at run time by `YAFU_NEON_MED=ranges|split|portable`.
  - N1, `tools/neon_bench.sh 'siqs(rsa(240))' 2` (C73, 1 thread, -noopt),
    QS seconds / Mcycles / polys, two reps each within 0.1 s:

    | config | QS s | Mcyc | polys |
    |---|---|---|---|
    | `-noneon` (all portable) | 48.2 | 123,115 | 113,929 |
    | ranges | 41.1 | 105,620 | 112,786 |
    | split | 39.6 | 102,040 | 114,564 |

    `-noneon` also turns off the NEON direct-sieve kernels and scanner, so
    the med sieve's own gain needs `YAFU_NEON_MED=portable` vs `split`.
- **P4 tdiv_med / resieve_med:** 16-bit Barrett checks with the
  `pmulhuw` replacement.
- **P5 tdiv_LP:** 4-lane compare of bucket entries.
- **P6** firstRoots, lp_sieveblock, SPV: only what the profile shows.
- **P7 Lanczos:** `v_t` as `uint64x2_t`.
- **P8 tuning:** a NEON column in the per-ISA params tables (siqs_tune);
  try 64 KB sieve blocks (the 64k code paths exist) since N1 has 64 KB L1d;
  revisit the mid-range strategy (bucket vs sparse) with NEON kernels in.
