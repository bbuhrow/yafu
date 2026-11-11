# arith/limb — one- and two-limb arithmetic

A small, fast, self-contained library for arithmetic on one 64-bit limb and
on two (128 bits), scalar and vectorized. Everything YAFU's tiny factoring and
PRP routines need below "arbitrary precision" lives here, in one copy.

It depends on **nothing else in YAFU** and on no external library (no GMP):
copy this directory into another project, add it to the include path, and
compile `limb1.c` and `limb2.c` if you use anything that is not inline.

| file | contents |
|---|---|
| `mp_platform.h` | compiler / OS / ISA detection, inline attribute macros, and one polyfill each for the wide-multiply, wide-divide, carry-chain, bit-scan and popcount primitives |
| `limb1.h` / `limb1.c` | single limb: modular add/sub, Montgomery multiply/square (several flavours), inverses, gcd, modexp, jacobi, sqrt mod p; 52-bit AVX-512 kernels |
| `limb2.h` / `limb2.c` | two limbs: plain add/sub/shift/multiply/square/divide, 128-bit Montgomery arithmetic and setup, 128-bit gcd; 104-bit AVX-512 kernels |
| `limb_selftest.c` | stand-alone self-test and micro-benchmark (no GMP, no YAFU) |

`limb2.h` includes `limb1.h`, which includes `mp_platform.h`.

## Conventions

- A limb is `uint64_t`. A two-limb value is a little-endian pair: `x[0]` low,
  `x[1]` high.
- Hot kernels are `static` force-inlined functions in the headers
  (`MP_FORCE_INLINE`). Larger or colder routines are ordinary functions in
  the `.c` files.
- Modular routines take **reduced** inputs (in `[0, n)`) and return reduced
  outputs unless the name says otherwise (`mulredc60` is lazy: `[0, 2n)`).
- Montgomery moduli must be odd. Two inverse conventions are in use, because
  both are fastest somewhere:
  - **positive** inverse `inv = n^-1 mod R` — `mulredc_pos`, `sqrredc_pos`,
    `mfma64`, and the 52-bit x8 kernels. Get it from
    `multiplicative_inverse(n)` / `multiplicative_inverse52(n)`.
  - **negative** inverse `nhat = -n^-1 mod R` — `mulredc`, `sqrredc`,
    `mulredc63`, `mulmod128` (`monty128_t.rho`). It is `0 - inv`.

## API index

### `mp_platform.h`

| | |
|---|---|
| 64 x 64 → 128 | `_umul128(a, b, &hi)`, `__umulh(a, b)` |
| 128 / 64 → 64 | `_udiv128(hi, lo, d, &rem)` (needs `hi < d`), `mp_udiv128_soft` (portable, same contract) |
| carry chains | `mp_addcarry_u64/_u32`, `mp_subborrow_u64/_u32` — return the carry, write the result through the pointer, on every compiler |
| bit scans | `_trail_zcnt64`, `_lead_zcnt64`, `_trail_zcnt` (undefined at 0); `_trail_full_zcnt`, `_lead_full_zcnt` (return 64 at 0); `_reset_lsb64`, `_reset_lsb` |
| counting | `mp_popcount64`, `mp_popcount32` |
| configuration | `MP_ARCH_X64`, `MP_ARCH_ARM64`, `MP_HAS_GNU_ASM_X64`, `MP_USE_BMI2`, `MP_USE_AVX2`, `MP_USE_AVX512F`, `MP_USE_IFMA`, `HAS_UINT128` |

### `limb1.h` (inline)

| | |
|---|---|
| modular add / sub | `addmod`, `submod`, `addmod32`, `submod32` |
| Montgomery, positive inverse | `mulredc_pos`, `sqrredc_pos`, `mulredc_pos_alt`, `mfma64` (multiply-add) |
| Montgomery, negative inverse | `mulredc`, `sqrredc` (any odd n), `mulredc63`, `sqrredc63` (n < 2^63), `mulredc60`, `sqrredc60` (n < 2^60, lazy) |
| inverses mod 2^k | `multiplicative_inverse` (2^64), `multiplicative_inverse52` (2^52) |
| gcd | `bingcd64(u, v)` (u odd) |
| 52-bit, 8 lanes (AVX-512) | `mulredc52_x8`, `sqrredc52_x8`, `addmod52_x8`, `submod52_x8`, `mp_mul52lo_x8`, `mp_mul52hi_x8`, `mp_mul52lohi_x8` |
| 52-bit vector building blocks | `mul52lo/hi/lohi`, `VEC_MUL_ACCUM_LOHI_PD`, `_mm512_{addsetc,subsetc,adc,sbb,addcarry,subborrow}_epi52` and masked forms, `carryprop` |

### `limb1.c`

`spAdd`, `spAdd3`, `spSub`, `spSub3`, `spMultiply`, `spMulAdd`, `spMulMod`,
`spDivide`, `u64div`, `spModExp`, `spPRP2`, `spBits`, `bits64`, `ndigits_1`,
`modinv_1`, `modinv_1b`, `modinv_1c`, `jacobi_1`, `ShanksTonelli_1`, `spGCD`,
`gcd64`, `spBinGCD`, `spBinGCD_odd`, `dblGCD`.

### `limb2.h` / `limb2.c`

| | |
|---|---|
| plain arithmetic (inline) | `add128`, `sub128`, `cmp128`, `shl128`, `shr128`, `mul128x64` (→192), `mul128` (→256), `sqr128` (→256) |
| division | `uint128_div(n, d, q, r)`, `mod_128(a, b, c)` |
| Montgomery setup | `monty128_init`, `to_monty128` |
| Montgomery arithmetic | `mulmod128`, `sqrmod128`, `mulmod128n`, `sqrmod128n`, `addmod128`, `submod128`, `dblmod128`, `chkmod128` |
| bit scans | `my_clz128`, `my_ctz128`, `my_clz104`, `my_ctz104` (zero-safe) |
| gcd | `gcd128` (all compilers), `bin_gcd128` (compilers with `__int128`; u odd, clobbers its inputs) |
| 104-bit, 8 lanes (AVX-512) | `mask_mulredc104_vec`, `mulredc104_vec`, `sqrredc104_vec`, `mask_sqrredc104_vec`, `mask_sqrredc104_vec_pos`, `mask_sqrredc104_exact_vec`, `addmod104_x8`, `submod104_x8`, `chkmod104_x8`, `mask_addmod104_x8`, `mask_dblmod104_x8`, `mask_redsub104_x8`, `redsub104_x8`, `multiplicative_inverse104_x8`, `mulredc52_mask_add_vec` |

## Build configuration

No flags are required. The library reads the compiler's own feature macros
(`__BMI2__`, `__AVX2__`, `__AVX512F__` + `__AVX512DQ__`, `__AVX512IFMA__`),
so `-march=native`, `-mavx2 -mbmi2`, `/arch:AVX512`, ... just work.

A host project can instead select features explicitly, the way YAFU does:
`-DUSE_BMI2`, `-DUSE_AVX2`, `-DUSE_AVX512F`, `-DIFMA`. **If any of those is
defined, auto-detection is switched off**, so the library always takes the
same path as the host's own `#ifdef`-gated code.

| macro | effect |
|---|---|
| `MP_NO_AUTO_ISA` | never infer features from compiler macros |
| `ASM_ARITH_DEBUG` | disable the GNU inline-asm paths (A/B against the portable code) |
| `MP_FORCE_PORTABLE_CARRY`, `MP_NO_NATIVE_ADDCARRY`, `MP_ASSUME_NATIVE_ADDCARRY` | choose the carry-chain implementation |
| `LIMB1_SCALAR_ONLY` | omit the 52-bit vector building blocks (short macro names such as `set64`, and file-scope statics) and the 104-bit layer; the `*_x8` kernels remain |
| `LIMB1_NO_VECTOR`, `LIMB2_NO_VECTOR` | omit all vector code of that header |

## Portability

Checked with `limb_selftest.c` (every routine against bit-at-a-time
references) on:

| compiler | target | modes |
|---|---|---|
| gcc 11 | x86-64 Linux | baseline, AVX2+BMI2, AVX-512, AVX-512+IFMA, asm disabled, portable carries |
| gcc 11 | aarch64 Linux (qemu) | baseline |
| MSVC 19.44 | x64 Windows | baseline, `/arch:AVX2`, `/arch:AVX512` |
| MSVC 19.44 | ARM64 Windows | compile only |
| clang-cl (VS 2022 LLVM) | x64 Windows | baseline, AVX-512, AVX-512+IFMA |
| Intel icx 2023.1 | x64 Windows | baseline, AVX-512, AVX-512+IFMA |
| Intel icl (classic) 2021.9 | x64 Windows | baseline, AVX2, AVX-512 |

Not yet run: clang on Linux/macOS (Apple silicon), AMD hardware, MinGW. The
code paths those take are covered by the gcc-x86, gcc-aarch64 and clang-cl
runs, but the binaries themselves have not been executed there.

What each platform gets:

| | x86-64 gcc / clang / icx / icl | x86-64 MSVC | ARM64 |
|---|---|---|---|
| 64-bit Montgomery, modular add/sub | GNU inline asm | intrinsics (`_umul128`, `_addcarry_u64`) | `__int128` (gcc/clang), `__umulh` (MSVC) |
| 128-bit Montgomery multiply | mulx asm (BMI2) or `__int128` CIOS | intrinsic CIOS | `__int128` CIOS |
| 128 / 64 divide | `divq` | `_udiv128` | `__int128` (gcc/clang), software (MSVC) |
| 52 / 104-bit vector kernels | AVX-512, IFMA when available | AVX-512 | — |

Vector kernels exist for AVX-512 only (8 lanes). There are no AVX2 (4-lane)
or NEON versions in this library; YAFU's AVX2 x4 ECM code keeps its own
double-precision arithmetic in `factor/gmp-ecm/edw_kern.h`.

## Testing

```
gcc -O2 -march=native limb_selftest.c limb1.c limb2.c -lm -o limb_selftest
./limb_selftest          # exit status 0 = all checks passed
./limb_selftest bench    # plus ns/op for the main kernels
```

Inside YAFU, `make test-run` also runs the `limb` module
(`test/layer1/test_limb.c`), which checks the same routines against GMP.
