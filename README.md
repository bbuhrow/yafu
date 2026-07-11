# YAFU — Yet Another Factoring Utility

Automated integer factorization. 

YAFU (with assistance from other free software) uses the most powerful modern algorithms (and implementations of them) to factor input integers in a completely automated way.  

YAFU has been referenced several times in the academic literature.  If you have academic work that requires integer factorization, YAFU might be able to help.  
 
For factorization help, support, and discussion, see the community at <https://www.mersenneforum.org/node/58>.

> **Detailed documentation lives in the [project wiki](../../wiki).** The wiki covers every option, every callable function, and a fuller version of the build instructions below.

---

## Quick start

### Use a pre-built Windows binary

Download from the [latest release](../../releases/latest). All Windows builds are **fully static** — no extra DLLs required, runs in `cmd.exe` and PowerShell.

| File | Best for |
| --- | --- |
| `yafu-windows-generic.exe`    | Any x86-64 CPU — widest compatibility |
| `yafu-windows-sse41.exe`      | Intel Core 2 / Penryn (2008) and newer |
| `yafu-windows-avx2.exe`       | Haswell (2013) and newer — **recommended for most users** |
| `yafu-windows-avx512.exe`     | Skylake-X / Cascade Lake / Zen 4 (AVX-512) |
| `yafu-windows-avx512ifma.exe` | Ice Lake and newer — best ECM performance |

Not sure which to pick? Run `wmic cpu get name` in a Windows command prompt and look up your CPU generation, or just start with `avx2` — it works on nearly every CPU sold since 2015.

### First run

```
yafu "factor(rsa(200))"
```

Or launch the interactive prompt:

```
yafu
>> factor(2056802480868100646375721251575555494408897387375737955882170045672576386016591560879707933101909539325829251496440620798637813)
```

Files YAFU writes (logs, savefiles, `yafu.ini`) go next to the executable. No installation step is needed.

---

## Building from source

Pre-built Windows binaries are provided, but **building yourself, tuned to your CPU, will very likely be faster.** The whole build is driven by a `Makefile`, configured per-machine through `config.mk`.

### Dependencies

| Required | Optional |
| --- | --- |
| **GMP** — <https://gmplib.org/> | **GMP-ECM** — <http://ecm.gforge.inria.fr/> (for ECM factorization) |
| | **CUDA Toolkit** (for GPU cofactorization / poly select / LA) |
| | **OpenCL** (alternative to CUDA for cofactorization) |

As of YAFU 3.0, `ytools`, `ysieve`, and `msieve` are bundled in this repository and built as part of the YAFU build.

### Linux, WSL, and MSYS2/MinGW-w64

```bash
cp config.mk.example config.mk          # first time only — edit paths as needed
make yafu                               # builds with sensible defaults (gcc)
```

Pass feature/ISA flags either on the command line or by uncommenting them in `config.mk`. The most useful ones are listed below. For the complete set, see the [Building YAFU](../../wiki/Building-YAFU) wiki page or run `make help`.

```bash
# Inspect the resolved configuration before/during a build
make info
make help

# Using config.mk for configuration of ISA and features
make yafu

or

make all

# ISA selection — overriding config.mk with specified ISA
make yafu USE_AVX2=1
make yafu USE_AVX512=1          # Skylake-X, Cascade Lake, Zen 4, etc.
make yafu USE_AVX512IFMA=1      # Ice Lake and newer

# ISA selection — overriding config.mk with feature flags or compiler options
make yafu ECM=1 OMP=1           # link GMP-ECM, enable OpenMP
make yafu CC=clang              # use clang
make yafu DEBUG=1               # debug build

```

The first time you set up, copy `config.mk.example` to `config.mk` and edit any paths that differ from your system. `config.mk` is gitignored — only `config.mk.example` is tracked.

### Other useful targets

| Target | Builds |
| --- | --- |
| `yafu` | The main yafu executable |
| `msieve` | The msieve static library + demo |
| `siqs` | The siqs_demo standalone binary |
| `ecm` | The ecm_demo standalone binary |
| `all` | All four above |
| `info` | Print the fully resolved configuration |
| `help` | Show the feature-flag reference |
| `clean` | Remove build artifacts |

### Windows — Visual Studio (MSVC)

Build files live under `build.vc22/`. To configure your build:

1. Open `build.vc22/yafu.sln` in Visual Studio 2022 (the free Community edition is fine).
2. Edit the property sheets in `build.vc22/` for your environment:
   - **`Directory.Build.props`** — global build settings inherited by every project.
   - **`build_config.props`** — toolset selection (e.g. v143, Intel oneAPI) and global build configuration.
   - **`libs_and_extensions.props`** — paths to the dependency libraries and the ISA extensions to enable for this build.
3. Build the `yafu` project under `x64 / Release`.

Dependencies for the MSVC build:

| Dependency | Source |
| --- | --- |
| **mpir** (or GMP) | <https://github.com/BrianGladman/mpir> |
| **ecm** | <https://gitlab.inria.fr/zimmerma/ecm/-/releases> or <https://github.com/sethtroisi/gmp-ecm> |
| **pthreads** | <https://github.com/BrianGladman/pthreads> |

These ship with their own MSVC build files. Build them first; then point YAFU at them via the paths in `libs_and_extensions.props`.

The free Intel compiler also works inside Visual Studio 2022.

---

## Continuous integration & release builds

YAFU's release artifacts (the matrix of Windows `.exe` files shown above) are produced by an **automated GitHub Actions pipeline** that exercises:

- **Linux / gcc**, via `Makefile`
- **Windows / MSYS2 / MinGW-w64 gcc**, via `Makefile` with `STATIC_WIN=1` for the fully-static binaries
- **Windows / MSVC**, via the `build.vc22/` solution and property sheets (this leg is being expanded for CI)

For every push, the CI builds the full matrix of ISA targets (`generic`, `sse41`, `avx2`, `avx512`, `avx512ifma`) on each platform. On tagged releases (e.g. `v3.1.2`), the resulting binaries are attached to the release page.

See the [Continuous Integration](../../wiki/Continuous-Integration) wiki page for details on the matrix, the artifacts, and how the same `Makefile` flags drive both CI and local builds.

---

## GGNFS sievers (required for NFS)

For NFS factorizations, YAFU needs external GGNFS lattice sieve binaries (`ggnfs-lasieve4I*`). Linux and MinGW binaries are bundled under `factor/lasieve5_64/bin/`. Point YAFU at them with `ggnfs_dir=` in `yafu.ini`, or `-ggnfs_dir <path>` on the command line. Without these, NFS will not run.

If you have AVX-512 on your CPU, YAFU will also use **AVX-ECM** as the default ECM backend (built in). A standalone version lives at <https://github.com/bbuhrow/avx-ecm>.

---

## Help and documentation

| Where | What |
| --- | --- |
| `help` at the YAFU prompt | Built-in help (reads `docfile.txt`) |
| `help <function>` | Per-function detail |
| `docfile.txt` | Function reference (must sit next to the binary for `help` to work) |
| `yafu.ini` | Every option, as commented lines |
| [Project wiki](../../wiki) | All of the above, cross-referenced and searchable |
| <https://www.mersenneforum.org/node/58> | Community / support |

---

## Fun examples

```
# Exercises many of yafu's algorithms
yafu "factor(2056802480868100646375721251575555494408897387375737955882170045672576386016591560879707933101909539325829251496440620798637813)"

# Neat ECM example
yafu "ecm(140870298550359924914704160737419905257747544866892632000062896476968602578482966342704)"
```

If you build YAFU on a new platform or with a new compiler — or build a binary that beats one of the pre-compiled releases for a specific CPU — the maintainer would like to hear about it.
