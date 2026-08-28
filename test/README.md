# YAFU modular test system

A single, layered test runner for YAFU, replacing the scattered ad-hoc
`test_*` routines with one harness that reports a pass/fail table and exits
non-zero on any failure (CI-wireable).

## Layout

```
test/
  testkit.h / testkit.c     the framework: registration, assertions,
                            deterministic RNG, micro-benchmarks, CLI runner
  test_data.h / test_data.c shared known-answer corpus (primes, composites,
                            strong pseudoprimes, Carmichaels, semiprimes) +
                            an exact uint64 primality oracle + generators,
                            plus a corpus self-check module
  test_main.c               lists the modules to run
  layer0/                   primitive tests (mp_platform.h)
    test_mp_arith.c           umul128 / udiv128 / add-/sub-carry chains
    test_mp_bitscan.c         lead/trail zero-count, reset-lsb
  layer1/                   kernel tests (link YAFU + GMP)
    test_sp_arith.c           arith.h single-word mul/div/mulmod/modexp/modinv
    test_modular.c            Montgomery (mpz + 128-bit) vs GMP
    test_primality.c          APRCL/BPSW/MR + tinyprp vs the corpus
  Makefile
```

## The layer model

- **Layer 0 — primitives.** The platform arithmetic in `mp_platform.h`.
  Checked against self-contained portable-C reference oracles; no external
  dependencies. This is where the `mp_arith` / `mp_bitscan` harnesses now live.
- **Layer 1 — kernels.** Arithmetic and number-theory building blocks
  (single-word arith, Montgomery, primality). GMP is the reference oracle;
  these link against YAFU objects.
- **Layer 2 — factorization routines** (rho, squfof, ecm, siqs, tinyqs, …)
  and **Layer 3 — integration** (the `factor()` dispatcher, `calc`, black-box
  binary) are the next passes and are not included here.

Adding a module is: write `t_*` test functions taking `tk_ctx *`, list them
in a `tk_test[]`, expose one `const tk_module`, and add its `extern` to
`test_main.c`. No registration macros or linker tricks.

## Building

Default (framework + Layer 0). Needs only your include path so it can find
`mp_platform.h`; no GMP, no YAFU objects:

```
make subset            # or just: make
./yafu_test
```

Full (adds Layer 1). Point it at your tree and the YAFU objects to link:

```
make full YAFU_ROOT=/path/to/yafu YAFU_OBJS="...objects except the driver..."
./yafu_test
```

See the Makefile header for the variables. Layer 1 is gated behind
`-DTK_NO_LAYER1` (the subset build defines it), so the framework and Layer 0
always build even without GMP.

## Running

```
./yafu_test                     run everything
./yafu_test --list              list modules and tests with their tags
./yafu_test mp_arith modular    run only these modules (bare args = names)
./yafu_test --tag primitives    run only tests carrying a tag
./yafu_test --seed 12345        set the RNG seed (printed so failures reproduce)
./yafu_test --bench             also run micro-benchmarks (never affects pass/fail)
./yafu_test --stop              stop after the first failing test
```

Assertions are non-aborting by default (`TK_CHECK*`): a test reports every
failure it finds in one run. `TK_REQUIRE` aborts just the current test when
continuing is pointless. Each test gets an independent, reproducible RNG
stream derived from the seed.

## Status of this pass (for review)

- **Framework, shared data, Layer 0** are compile-clean (`-Wall -Wextra
  -Wpedantic`) and were run: the corpus self-check, `meta`, `mp_arith`, and
  `mp_bitscan` all pass. Layer 0 was exercised here against a portable stand-in
  for `mp_platform.h`; in your tree it includes the real header unchanged, so
  it then measures the actual intrinsics (visible under `--bench`).
- **Layer 1** is written against the real YAFU/GMP signatures and passes a
  `-Wall -Wextra` syntax check against those prototypes, but has **not** been
  run here (no GMP/YAFU objects in the scratch environment). Expect to build
  it in-tree first. A few call conventions I couldn't compile against are
  flagged inline with `NOTE:` / `[LIMB-ORDER]` / `[PRP-CONV]` — grep for those
  and confirm on first build.

### Deliberately deferred
`test.sh`, `testbench.txt`, all NFS, and the `lasieve*` sievers are out of
scope for now. Layer 2 (factorization routines, incl. the `microecm` harness)
and Layer 3 (integration) come after this is reviewed. A `sieve` module
(ysieve `soe` vs known π(x)) is the natural next Layer-1 addition.
