# YAFU GUI — prototype

A calculator-style desktop front end for YAFU, in Python + Qt (PySide6).

```
pip install PySide6
python yafu_gui.py --yafu /path/to/yafu --cwd /path/to/yafu_dir
```

`--yafu` and `--cwd` are optional. The GUI looks for yafu next to itself, in the parent folder and on PATH, and remembers your choice (File ▸ Preferences). The working directory is where yafu reads `yafu.ini` and writes `siqs.dat`, `nfs.dat` and `factor.log`. It defaults to the executable's folder.

## Files

| File | Purpose |
|---|---|
| `yafu_gui.py` | The application |
| `yafu_options.json` | All 149 options (name, help, needsArg, inferred type), generated from `cmdOptions.c` |
| `yafu_schema.py` | Regenerates the JSON: `python yafu_schema.py top/cmdParser/cmdOptions.c` |
| `yafu_ini.py` | Reads and writes `yafu.ini` with `readINI()` semantics. Only the changed lines are touched, and comments are kept |
| `driver_yafu_gui.patch` | Optional 2-hunk `driver.c` patch (see "Platform notes") |
| `tests/fake_yafu.py` | A stand-in yafu that copies driver.c's I/O behaviour, for testing |
| `tests/test_gui.py` | Headless end-to-end test (25 checks): `QT_QPA_PLATFORM=offscreen python tests/test_gui.py` |

## How it works

- **YAFU is unchanged.** Each job runs as one `yafu "<expr>" -opt val …` process, the same command you would type. The read-only line above the Run button shows that exact command, and you can copy it.
- **Modes**: factor(), SIQS, NFS, ECM, P±1, Rho/Fermat/trial, Primes (SoE), Primality and Calculator. Each mode wraps the input in its function (`siqs(x)`, `ecm(x,curves)`, `pp1(x,bases)`, `primes(lo,hi[,0])`, `aprcl(x)` …). The mode decides which option groups are expanded, which are collapsed, and which are hidden.
  - factor() shows every algorithm group (collapsed), because factor() uses all of them.
  - Options in hidden groups are not passed. The tooltip on the command line lists any that are set but not applied in the current mode.
- **Defaults layering matches yafu**: built-in `initOpt()` defaults < `yafu.ini` < command line.
  - Each field's placeholder shows the effective default and where it comes from, e.g. `11000 (yafu.ini)` or `auto`.
  - Only options you change are passed. Changed options are highlighted, and each group header counts its changes.
- **Settings ▸ Defaults (yafu.ini)** edits the ini as a table of every option. Saving writes a `.bak` copy and changes only the affected lines. Disabled options are commented out; enabled ones uncomment the existing example line where there is one.
- **Multi-line input**: in an algorithm mode it runs as `-batchfile` with `func(@)`. In Calculator mode it runs as `-script`.
- **Results pane**: parsed from the console output (`Pnn = …`, `***factorization:***`, `ans = …`) and grouped per input. Batches split on yafu's `=== Starting work on batchfile expression #N ===` marker.
- **Other features**:
  - Size preview (`C155 · 514 bits`) for plain arithmetic input.
  - ECM B1 presets (t20–t65).
  - History tab: double-click to restore a job.
  - Ctrl+Enter runs the job; Ctrl+↑/↓ recall earlier inputs.
- **Stop** sends Ctrl-C (SIGINT) to the whole process group. SIQS's handler saves `siqs.dat`, and external sievers are interrupted too. A second click kills the process.
- **Future-proofing**: options added to `cmdOptions.c` appear automatically under "Other options" after re-running `yafu_schema.py`. The curated layout in `MODE_GROUPS` only refers to option names.

## Platform notes (from reading `driver.c`)

Two `driver.c` behaviours break a naive subprocess wrapper:

1. `check_expression()` treats a non-TTY stdin as **"batchfile from stdin"**. If stdin is a pipe or `/dev/null`, yafu ignores the command-line expression.
2. When stdout is a pipe, the C runtime **fully buffers** it, so `\r` progress lines arrive in 4 KB bursts or only at exit.

- **Linux/macOS**: the GUI runs yafu on a pseudo-terminal, so `isatty()` is true and stdout is line-buffered. Stock yafu works with no changes.
- **Windows**: the GUI uses pipes, with stdin connected to `NUL`. MSVC's `_isatty(NUL)` returns true, and mingw builds skip pipe detection, so jobs run correctly but console output arrives in chunks. `driver_yafu_gui.patch` fixes this when the `YAFU_GUI` environment variable is set (the GUI always sets it). It makes stdout unbuffered and skips pipe detection, and has no effect otherwise. With the patch, "pipes" also works on POSIX (Preferences ▸ Process I/O).
- **Windows Stop**: there is no portable way to send Ctrl-C to a child console process, so Stop kills it. A possible follow-up: start yafu with `CREATE_NEW_PROCESS_GROUP` and send `CTRL_BREAK_EVENT`, which would need a SIGBREAK handler beside the SIGINT one in SIQS.

## Limitations

- **Flags set in `yafu.ini` can't be turned off from the command line**, because yafu has no `-no<flag>`. The GUI shows them checked and locked, with a tooltip pointing to the Defaults editor. A `-ini <file>` option in yafu would remove this limit: the GUI could write a complete per-job ini.
- **Verbosity** can only be raised above the ini level.
- **Not tested against the real binary here.** The repo's `bin/linux/yafu` is an icc build that needs `libsvml.so`. Everything was tested against `fake_yafu.py`, which copies the relevant driver.c behaviour.

## Found while reading the source

- **Quotes in `yafu.ini` values are kept literally.** `readINI()` keeps quotes in values, but the shipped `yafu.ini` comments show examples like `logfile="factor.log"` and `script="my_script_file"`. Uncommenting one gives a filename that includes the quote characters, which fails outright on Windows. The GUI writes values without quotes.
- **`-B1ecm` is the only numeric option that accepts `1e6`-style values.** Every other numeric option goes through `enforce_numeric()`, which rejects that notation. The GUI's validators follow this rule.
