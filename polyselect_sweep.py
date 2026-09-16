#!/usr/bin/env python3
"""
polyselect_sweep.py -- run YAFU NFS poly select over a grid of
(input, coeff_deadline, repeat), capture each run's stdout+stderr to a log,
parse the logs, and tabulate the results.

Parsed per run:
  hits     last "| hits X" value from the status roll-up (k/M/G expanded)
  cand     last "cand X"
  saved    last "saved X"
  ad_idx   last "a_d N (i/T)" -> i and T
  spq      last "spq: X/Y" -> X
  elapsed  "elapsed time: S seconds"
  best_e   "e <value>" on the "# norm ..." line(s) (max if several)

Examples
  # run the sweep (serial; each run gets a fresh working dir)
  python3 polyselect_sweep.py run --yafu /path/to/yafu --ini /path/to/yafu.ini \
      --inputs inputs.txt --deadlines 3,5,15,30 --threads 16

  # re-parse existing logs only
  python3 polyselect_sweep.py parse --outdir sweep_out

  # parse arbitrary log files
  python3 polyselect_sweep.py parse capture.txt
"""

import argparse
import csv
import os
import pty
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from collections import defaultdict

# --------------------------------------------------------------------------
# Command construction -- ADJUST to match how your build takes the args.
#   {yafu} {threads} {cd} {n} {rep} are substituted.
# The expression goes on the command line; stdin is a pty (see run_sweep).
# --------------------------------------------------------------------------
DEFAULT_CMD = ('{yafu} "nfs({n})" -np -threads {threads} -v '
               '-nfs_stage1_args "coeff_deadline={cd}"')

SUFFIX = {"": 1, "k": 1e3, "K": 1e3, "M": 1e6, "G": 1e9, "T": 1e12}
NUM = r"([0-9]*\.?[0-9]+)\s*([kKMGT]?)"
FLT = r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"

RE_HITS = re.compile(r"\|\s*hits\s+" + NUM)
RE_CAND = re.compile(r"\bcand\s+" + NUM)
RE_SAVED = re.compile(r"\bsaved\s+" + NUM)
RE_AD = re.compile(r"a_d\s+(\d+)\s+\((\d+)/(\d+)\)")
RE_SPQ = re.compile(r"spq:\s*" + NUM + r"\s*/\s*" + NUM)
RE_ELAPSED = re.compile(r"^elapsed time:\s*" + FLT +
                        r"(?:\s*seconds\s*\(\s*" + FLT + r"\s*second deadline)?",
                        re.M)
RE_NORMLINE = re.compile(r"^#\s*norm\b.*$", re.M)
RE_E = re.compile(r"\be\s+" + FLT)
RE_BEST_AD = re.compile(r"best E\s+" + FLT + r"\s*@\s*a_d\s+(\d+)")


EARLY_MARGIN = 1.0   # seconds; set via --early-margin


def num(m, i=1):
    return float(m.group(i)) * SUFFIX[m.group(i + 1)]


def last(regex, text):
    ms = list(regex.finditer(text))
    return ms[-1] if ms else None


def parse_log(text):
    """Return dict of parsed values (None where not found)."""
    r = dict(hits=None, cand=None, saved=None, ad=None, ad_idx=None,
             ad_total=None, spq=None, spq_total=None, elapsed=None,
             best_e=None, best_e_ad=None, deadline=None, early=None)
    # status roll-up may be \r-separated; regexes work on the raw text
    m = last(RE_HITS, text)
    if m: r["hits"] = int(num(m))
    m = last(RE_CAND, text)
    if m: r["cand"] = int(num(m))
    m = last(RE_SAVED, text)
    if m: r["saved"] = int(num(m))
    m = last(RE_AD, text)
    if m:
        r["ad"], r["ad_idx"], r["ad_total"] = (int(g) for g in m.groups())
    m = last(RE_SPQ, text)
    if m:
        r["spq"] = int(num(m, 1))
        r["spq_total"] = int(num(m, 3))
    m = last(RE_ELAPSED, text)
    if m:
        r["elapsed"] = float(m.group(1))
        if m.group(2):
            r["deadline"] = float(m.group(2))
            # early exit = finished comfortably before the poly deadline
            r["early"] = int(r["elapsed"] < r["deadline"] - EARLY_MARGIN)
    es = []
    for line in RE_NORMLINE.finditer(text):
        me = RE_E.search(line.group(0))
        if me:
            es.append(float(me.group(1)))
    if es:
        r["best_e"] = max(es)
    m = last(RE_BEST_AD, text)
    if m: r["best_e_ad"] = int(m.group(2))
    return r


# --------------------------------------------------------------------------
def read_inputs(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                out.append(line)
    return out


def log_name(idx, n, cd, rep):
    return f"in{idx:03d}_c{len(n)}_cd{cd}_r{rep}.log"


def run_sweep(a):
    inputs = read_inputs(a.inputs)
    os.makedirs(a.outdir, exist_ok=True)
    logdir = os.path.join(a.outdir, "logs")
    os.makedirs(logdir, exist_ok=True)
    yafu = os.path.abspath(a.yafu)
    # accept "3 5 15 30", "3,5,15,30", or a mix
    a.deadlines = [d for tok in a.deadlines for d in tok.split(",") if d]
    ini = os.path.abspath(a.ini) if a.ini else None

    jobs = [(i, n, cd, rep) for i, n in enumerate(inputs)
            for cd in a.deadlines for rep in range(a.repeats)]
    print(f"{len(jobs)} runs -> {a.outdir}", file=sys.stderr)

    for k, (i, n, cd, rep) in enumerate(jobs, 1):
        logpath = os.path.join(logdir, log_name(i, n, cd, rep))
        if a.resume and os.path.exists(logpath):
            print(f"[{k}/{len(jobs)}] skip {os.path.basename(logpath)}",
                  file=sys.stderr)
            continue

        # fresh working dir so no stale nfs.dat / .p / .job gets reused
        wd = os.path.join(a.outdir, "work", f"in{i:03d}_cd{cd}_r{rep}")
        shutil.rmtree(wd, ignore_errors=True)
        os.makedirs(wd)
        if ini:
            shutil.copy(ini, os.path.join(wd, "yafu.ini"))

        fmt = dict(yafu=shlex.quote(yafu), threads=a.threads, cd=cd,
                   n=n, rep=rep)
        # run the template through the shell verbatim, exactly as typed
        cmd = a.cmd.format(**fmt)
        shown = cmd

        print(f"[{k}/{len(jobs)}] in{i} c{len(n)} cd={cd} rep={rep}: "
              f"{shown}", file=sys.stderr)
        t0 = time.time()
        with open(logpath + ".tmp", "w") as lf:
            lf.write(f"### cmd: {shown}\n")
            lf.flush()
            # yafu treats a non-tty stdin as a batchfile, so give it a
            # pseudo-terminal (kept open until the run finishes)
            master, slave = pty.openpty()
            try:
                subprocess.run(cmd, shell=True, executable="/bin/bash",
                               cwd=wd, stdin=slave,
                               stdout=lf, stderr=subprocess.STDOUT,
                               timeout=a.timeout)
            except subprocess.TimeoutExpired:
                lf.write("\n### TIMEOUT\n")
            finally:
                os.close(slave)
                os.close(master)
        os.replace(logpath + ".tmp", logpath)
        print(f"    wall {time.time() - t0:.1f}s", file=sys.stderr)
        if not a.keep_work:
            shutil.rmtree(wd, ignore_errors=True)

    return sorted(os.path.join(logdir, f) for f in os.listdir(logdir)
                  if f.endswith(".log"))


# --------------------------------------------------------------------------
RE_NAME = re.compile(r"in(\d+)_c(\d+)_cd([^_]+)_r(\d+)\.log$")
COLS = ["input", "digits", "cd", "rep", "hits", "cand", "saved",
        "ad_idx", "ad_total", "spq", "elapsed", "deadline", "early", "best_e",
        "best_e_ad", "log"]


def collect(paths):
    rows = []
    for p in paths:
        with open(p, errors="replace") as f:
            r = parse_log(f.read())
        m = RE_NAME.search(os.path.basename(p))
        if m:
            r.update(input=int(m.group(1)), digits=int(m.group(2)),
                     cd=m.group(3), rep=int(m.group(4)))
        else:
            r.update(input="", digits="", cd="", rep="")
        r["log"] = os.path.basename(p)
        rows.append(r)
    return rows


def fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.4g}" if abs(v) < 1e-3 else f"{v:.2f}"
    return str(v)


def print_table(rows, cols, out=sys.stdout):
    cells = [[fmt(r.get(c)) for c in cols] for r in rows]
    w = [max(len(c), *(len(x[i]) for x in cells)) if cells else len(c)
         for i, c in enumerate(cols)]
    print("  ".join(c.rjust(w[i]) for i, c in enumerate(cols)), file=out)
    print("  ".join("-" * w[i] for i in range(len(cols))), file=out)
    for x in cells:
        print("  ".join(v.rjust(w[i]) for i, v in enumerate(x)), file=out)


def cd_key(v):
    try:
        return (0, float(v))
    except (TypeError, ValueError):
        return (1, str(v))


def summarize(rows):
    """Group by (input, cd): mean/median/max best E, means of the rest."""
    g = defaultdict(list)
    for r in rows:
        g[(r["input"], r["digits"], r["cd"])].append(r)
    out = []
    for (inp, dig, cd), rs in sorted(g.items(),
                                     key=lambda kv: (str(kv[0][0]),
                                                     cd_key(kv[0][2]))):
        def vals(k):
            return [r[k] for r in rs if r.get(k) is not None]

        def mean(k):
            v = vals(k)
            return statistics.fmean(v) if v else None
        e = vals("best_e")
        ea = [r for r in rs if r.get("early") is not None]
        early = [r["elapsed"] for r in ea if r["early"]]
        out.append(dict(
            input=inp, digits=dig, cd=cd, n=len(rs),
            e_mean=statistics.fmean(e) if e else None,
            e_med=statistics.median(e) if e else None,
            e_max=max(e) if e else None,
            e_sd=statistics.stdev(e) if len(e) > 1 else None,
            hits=mean("hits"), saved=mean("saved"),
            ad_idx=mean("ad_idx"), elapsed=mean("elapsed"),
            early_frac=(len(early) / len(ea)) if ea else None,
            early_t=statistics.fmean(early) if early else None))
    return out


SUMCOLS = ["input", "digits", "cd", "n", "e_mean", "e_med", "e_max", "e_sd",
           "hits", "saved", "ad_idx", "elapsed", "early_frac", "early_t"]


def report(rows, outdir):
    rows.sort(key=lambda r: (str(r["input"]), cd_key(r["cd"]),
                             str(r["rep"])))
    print("\n== per run ==")
    print_table(rows, COLS)
    summ = summarize(rows)
    print("\n== summary by input, coeff_deadline ==")
    print_table(summ, SUMCOLS)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        for name, data, cols in (("runs.csv", rows, COLS),
                                 ("summary.csv", summ, SUMCOLS)):
            with open(os.path.join(outdir, name), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
                w.writeheader()
                w.writerows(data)
        print(f"\nwrote {outdir}/runs.csv, {outdir}/summary.csv",
              file=sys.stderr)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    r = sub.add_parser("run", help="run the sweep, then parse")
    r.add_argument("--yafu", required=True, help="path to yafu binary")
    r.add_argument("--ini", help="yafu.ini to copy into each work dir")
    r.add_argument("--inputs", required=True,
                   help="file with one input (number/expression) per line")
    r.add_argument("--deadlines", nargs="+", required=True,
                   help="coeff_deadline values (space- or comma-separated)")
    r.add_argument("--repeats", type=int, default=1)
    r.add_argument("--threads", type=int, default=os.cpu_count())
    r.add_argument("--cmd", default=DEFAULT_CMD,
                   help=f"command template (default: {DEFAULT_CMD!r})")
    r.add_argument("--timeout", type=float, default=None,
                   help="per-run timeout, seconds")
    r.add_argument("--outdir", default="sweep_out")
    r.add_argument("--resume", action="store_true",
                   help="skip runs whose log already exists")
    r.add_argument("--keep-work", action="store_true",
                   help="keep per-run work dirs (nfs.dat.*.p etc.)")

    p = sub.add_parser("parse", help="parse existing logs only")
    p.add_argument("logs", nargs="*", help="log files (default: OUTDIR/logs/*.log)")
    p.add_argument("--outdir", default="sweep_out")

    for sp in (r, p):
        sp.add_argument("--early-margin", type=float, default=1.0,
                        help="run counts as early exit if elapsed < "
                             "deadline - margin (default 1.0 s)")

    a = ap.parse_args()
    global EARLY_MARGIN
    EARLY_MARGIN = a.early_margin
    if a.mode == "run":
        paths = run_sweep(a)
    else:
        paths = a.logs or sorted(
            os.path.join(a.outdir, "logs", f)
            for f in os.listdir(os.path.join(a.outdir, "logs"))
            if f.endswith(".log"))
    report(collect(paths), a.outdir)


if __name__ == "__main__":
    main()
