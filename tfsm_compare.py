#!/usr/bin/env python3
"""tfsm_compare.py -- is yafu's adaptive TFSm tuning worth it?

Compares, on the same inputs and binary with everything else at yafu's
defaults: no flags (adaptive TFSm), -noopt (the formula TFSm, no tuning), and
TFSm fixed at formula -8, -6, -4, -2, +2, +4.  4 balanced semiprimes of the
given size x 3 interleaved rounds (the config order rotates each round).
Runs are cached in <outdir>/runs.jsonl, so an interrupted or extended run
only does the missing work.

usage: tfsm_compare.py run|analyze <outdir> <digits> <threads>
  e.g. python3 tfsm_compare.py run /tmp/tfsm90 90 4 > /tmp/tfsm90.log 2>&1
       python3 tfsm_compare.py analyze /tmp/tfsm90 90 4

Uses the yafu binary and siqs_tune.py next to this script; Linux/WSL only
(each run gets a pseudo-terminal).
"""
import json
import math
import os
import pty
import re
import shutil
import statistics
import subprocess
import sys

REPO = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO)
import siqs_tune as st  # noqa: E402

YAFU = os.path.join(REPO, "yafu")
DIGITS = int(sys.argv[3]) if len(sys.argv) > 3 else 70
THREADS = int(sys.argv[4]) if len(sys.argv) > 4 else 1
NUMS = st.test_numbers(DIGITS, 4, "tfsm", "exp")
REPS = 3
DELTAS = [-8, -6, -4, -2, 2, 4]
RE_CHOSE = re.compile(r"SPV tuner: chose tf_small_cutoff = (\d+)")


def load(out):
    recs = []
    p = os.path.join(out, "runs.jsonl")
    if os.path.exists(p):
        for line in open(p):
            recs.append(json.loads(line))
    return recs


def run_one(out, n, cfg, extra, rep, done):
    key = f"{n}|{cfg}|{rep}"
    if key in done:
        return done[key]
    wd = os.path.join(out, "work")
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    cmd = [YAFU, f"siqs({n})", "-v", "-threads", str(THREADS), "-inmem", "1000"] + extra
    master, slave = pty.openpty()
    try:
        p = subprocess.run(cmd, cwd=wd, stdin=slave, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=600)
    finally:
        os.close(slave)
        os.close(master)
    text = p.stdout.decode(errors="replace")
    rec = st.parse_output(text, n)
    m = RE_CHOSE.search(text)
    rec.update(key=key, n=n, cfg=cfg, rep=rep, extra=extra,
               chosen=int(m.group(1)) if m else None)
    done[key] = rec
    with open(os.path.join(out, "runs.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"[{len(done)}] {n[:10]}.. {cfg:8s} rep{rep}: t={rec['t_total']} "
          f"TFSm={rec['TFSm']} chosen={rec['chosen']} partials={rec['partials']} "
          f"ok={rec['ok']}", flush=True)
    return rec


def cmd_run(out):
    os.makedirs(out, exist_ok=True)
    done = {r["key"]: r for r in load(out)}
    formula = {}
    for rep in range(REPS):
        for i, n in enumerate(NUMS):
            # -noopt first in round 0: its SPV correction is the formula value
            r = run_one(out, n, "noopt", ["-noopt"], rep, done)
            formula[n] = r["TFSm"]
            cfgs = [("adaptive", [])] + \
                [(f"F{d:+d}", ["-siqsTFSm", str(formula[n] + d)]) for d in DELTAS]
            k = (rep + i) % len(cfgs)          # rotate the order
            for cfg, extra in cfgs[k:] + cfgs[:k]:
                run_one(out, n, cfg, extra, rep, done)


def cmd_analyze(out):
    recs = [r for r in load(out)]
    bad = [r for r in recs if not r["ok"]]
    if bad:
        print(f"{len(bad)} failed runs: " + ", ".join(r["key"] for r in bad))
    recs = [r for r in recs if r["ok"]]
    cfgs = ["adaptive", "noopt"] + [f"F{d:+d}" for d in DELTAS]
    by = {}
    for r in recs:
        by.setdefault((r["cfg"], r["n"]), []).append(r)

    print(f"\nc{DIGITS}, {len(NUMS)} numbers x {REPS} repeats, {THREADS} thread(s)\n")
    print("per number: median total time [min-max], TFSm, median partials (spread)")
    for n in NUMS:
        f = by.get(("noopt", n), [{}])[0].get("TFSm")
        print(f"\n{n}  (formula TFSm = {f})")
        for c in cfgs:
            rs = by.get((c, n), [])
            if not rs:
                continue
            ts = [r["t_total"] for r in rs]
            ps = [r["partials"] for r in rs if r["partials"]]
            ch = [r["chosen"] for r in rs if r["chosen"] is not None]
            tf = sorted({r["TFSm"] for r in rs})
            pspread = (max(ps) - min(ps)) / statistics.mean(ps) if len(ps) > 1 else 0
            print(f"  {c:9s} t={statistics.median(ts):6.2f} [{min(ts):5.2f}-{max(ts):5.2f}]"
                  f"  TFSm={','.join(map(str, tf)) if not ch else 'chose ' + ','.join(map(str, ch))}"
                  f"  partials={statistics.median(ps):8.0f} ({100 * pspread:4.1f}%)")

    print("\nsummary vs adaptive (geometric mean over numbers of median-time ratio):")
    for c in cfgs:
        ratios = []
        for n in NUMS:
            a = by.get(("adaptive", n))
            b = by.get((c, n))
            if a and b:
                ratios.append(statistics.median(r["t_total"] for r in b) /
                              statistics.median(r["t_total"] for r in a))
        if ratios:
            g = math.exp(sum(map(math.log, ratios)) / len(ratios))
            # per-run noise: mean within-number relative spread of times
            spread = statistics.mean(
                (max(x["t_total"] for x in by[(c, n)]) - min(x["t_total"] for x in by[(c, n)])) /
                statistics.mean(x["t_total"] for x in by[(c, n)])
                for n in NUMS if (c, n) in by)
            print(f"  {c:9s} {g:6.3f}   per-number: " +
                  " ".join(f"{x:.3f}" for x in ratios) +
                  f"   (run-to-run time spread {100 * spread:.1f}%)")


if __name__ == "__main__":
    {"run": cmd_run, "analyze": cmd_analyze}[sys.argv[1]](sys.argv[2])
