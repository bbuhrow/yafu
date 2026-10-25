#!/usr/bin/env python3
"""
siqs_tune_v2.py -- two-pass SIQS parameter optimizer for yafu.

Pass 1 is a deterministic nested grid over the five parameters that matter
most, defined relative to what yafu itself chooses for each input, so the
results compare across inputs and sizes:

  B     factor base size, as a multiple of yafu's default B   0.5 0.65 0.8 1.0
  M     large prime multiplier                                 100 180 300
  TF    trial factoring bound, offset from yafu's default TF   -8 -4 0
  TFSm  small prime variation cutoff, offset from the formula
        value F (what yafu uses with -noopt)
                                     <= c70: -6 -4 -2   <= c80: -4 -2 0   else: -2 0 +2
  NB    sieve blocks (absolute)      <= c70: 1 2 3      <= c80: 2 3 4     else: 3 4 5

that is 4 x 3 x 3 x 3 x 3 = 324 points per input.  MFBD is fixed (1.72) and
so is the sieve pipeline (Mid/Large/Kernel, given explicitly so the NB sweep
can't switch pipelines underneath).  Pass 2 refines around the result.

Per digit size the workflow is:

  grid --inputs 1         every point on the first of 3 test numbers.  Each
                          number first gets a reference run (yafu's defaults
                          with -noopt and the fixed pipeline: gives the default
                          B and TF, the formula TFSm F and the reference time).
                          Points run in shuffled order; a run slower than
                          --timeout-mult x the reference is stopped and kept
                          as "slower than that".
  analyze                 per-level and pairwise tables, the best points, and
                          a suggested prune list (prune_cNN.txt) to review/edit
  grid --inputs 2,3 --prune
                          the surviving points on the other two numbers
  fit                     quadratic response surface in log B, log M, TF,
                          TFSm, NB over all numbers (per-number intercepts),
                          fitted near the best points: optimum, sensitivity
  trends                  the fitted optima across sizes, with linear trends
                          and predictions for larger sizes
  refine                  pass 2: small random steps around the fitted optimum
                          (MFBD and, at >= c70, Mid/Large/Kernel free again),
                          then a check on held-back numbers against yafu's
                          defaults and its built-in parameters on the bucket
                          pipeline -- written in siqs_tune.py's summary format
  export                  siqs_tune.py's export: siqs_candidates.txt for yafu's
                          siqstune()

Everything is cached in <outdir>/runs.jsonl (siqs_tune.py's format, so its
report/export work on a v2 outdir).  The cache doesn't know which yafu build
made a run, so <outdir>/binary.json records the binary and every command that
runs yafu refuses a different one (use a fresh outdir per build).

Examples
  python3 siqs_tune_v2.py grid --yafu ./yafu --digits 70 --inputs 1 --outdir v2out
  python3 siqs_tune_v2.py analyze --outdir v2out --digits 70
  python3 siqs_tune_v2.py grid --yafu ./yafu --digits 70 --inputs 2,3 --prune --outdir v2out
  python3 siqs_tune_v2.py fit --outdir v2out --digits 70
  python3 siqs_tune_v2.py trends --outdir v2out
  python3 siqs_tune_v2.py refine --yafu ./yafu --digits 70 --outdir v2out
  python3 siqs_tune_v2.py export --outdir v2out

Linux/WSL: yafu gets a pseudo-terminal per run (see siqs_tune.py).
"""

import argparse
import itertools
import json
import math
import os
import random
import statistics
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import siqs_tune as v1  # noqa: E402

GRID_PARAMS = ["B", "M", "TF", "TFSm", "NB"]
DEFAULT_PIPELINE = "Mid=sparse,Large=direct,Kernel=cmem"
FIT_TERMS = 5               # b, m, t, s, nb


def default_levels(d):
    if d <= 70:
        nb, ts = [1, 2, 3], [-6, -4, -2]
    elif d <= 80:
        nb, ts = [2, 3, 4], [-4, -2, 0]
    else:
        nb, ts = [3, 4, 5], [-2, 0, 2]
    return {"B": [0.5, 0.65, 0.8, 1.0], "M": [100, 180, 300],
            "TF": [-8, -4, 0], "TFSm": ts, "NB": nb}


def log(*args):
    print(*args, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------
# grid points: relative coordinates <-> absolute yafu parameters
# --------------------------------------------------------------------------
def fmt_level(k, v):
    if k == "B":
        return f"{v:g}"
    if k in ("TF", "TFSm"):
        return f"{int(v):+d}"
    return str(int(v))


def point_key(pt):
    return ",".join(f"{k}={fmt_level(k, pt[k])}" for k in GRID_PARAMS)


def parse_point(s):
    pt = {}
    for tok in s.replace(" ", ",").split(","):
        if not tok:
            continue
        k, v = tok.split("=", 1)
        if k not in GRID_PARAMS:
            raise ValueError(f"unknown grid parameter {k}")
        pt[k] = float(v) if k == "B" else int(v)
    return pt


def grid_points(levels):
    return [dict(zip(GRID_PARAMS, combo))
            for combo in itertools.product(*(levels[k] for k in GRID_PARAMS))]


def abs_params(pt, ref, meta):
    """yafu parameters for grid point pt on the number whose reference is ref"""
    p = dict(meta["pipeline"])
    p.update(B=v1.quantize("B", pt["B"] * ref["B"]),
             M=int(pt["M"]),
             TF=max(20, ref["TF"] + int(pt["TF"])),
             TFSm=max(1, ref["TFSm"] + int(pt["TFSm"])),
             NB=int(pt["NB"]),
             MFBD=meta["mfbd"])
    return p


def ref_params(meta):
    p = dict(meta["pipeline"])
    p["noopt"] = 1
    return p


# --------------------------------------------------------------------------
# per-size metadata and the binary guard
# --------------------------------------------------------------------------
def meta_path(outdir, d):
    return os.path.join(outdir, f"grid_c{d}.json")


def load_meta(outdir, d, required=True):
    p = meta_path(outdir, d)
    if not os.path.exists(p):
        if required:
            sys.exit(f"no {p}: run 'grid' for c{d} first")
        return None
    with open(p) as f:
        return json.load(f)


def save_meta(outdir, meta):
    with open(meta_path(outdir, meta["digits"]), "w") as f:
        json.dump(meta, f, indent=1)


def check_binary(a):
    yafu = os.path.abspath(a.yafu)
    st = os.stat(yafu)
    info = dict(path=yafu, size=st.st_size, mtime=int(st.st_mtime))
    p = os.path.join(a.outdir, "binary.json")
    os.makedirs(a.outdir, exist_ok=True)
    if os.path.exists(p):
        with open(p) as f:
            old = json.load(f)
        if (old["size"], old["mtime"]) != (info["size"], info["mtime"]):
            if not a.allow_binary_change:
                sys.exit(f"{yafu} is not the build this outdir was made with "
                         f"({old}); cached runs would mix builds.  Use a fresh "
                         f"--outdir, or --allow-binary-change if you mean it.")
            log(f"warning: binary changed since {old}")
    with open(p, "w") as f:
        json.dump(info, f, indent=1)


def make_runner(a):
    ns = argparse.Namespace(yafu=a.yafu, ini=a.ini, outdir=a.outdir,
                            keep_logs=a.keep_logs, threads=a.threads,
                            extra=a.extra, jobs=1)
    return v1.Runner(ns)


def init_meta(a, d):
    meta = load_meta(a.outdir, d, required=False)
    pipeline = {k: v1.quantize(k, v) for k, v in v1.parse_kv(a.pipeline).items()}
    if meta is None:
        meta = dict(digits=d, threads=a.threads, extra=a.extra,
                    pipeline=pipeline, mfbd=a.mfbd, seed=a.seed,
                    numbers=v1.test_numbers(d, a.numbers, a.seed, "grid"),
                    levels={k: [] for k in GRID_PARAMS}, refs={})
    else:
        for k, want in (("threads", a.threads), ("extra", a.extra),
                        ("pipeline", pipeline), ("mfbd", a.mfbd)):
            if meta[k] != want:
                sys.exit(f"c{d}: --{k} {want} differs from the existing grid's "
                         f"{meta[k]} in {meta_path(a.outdir, d)}")
    return meta


def reference(R, meta, n, reps, max_time):
    """reference runs for number n; cached in meta['refs']"""
    recs = [R.run(n, ref_params(meta), r, max_time) for r in range(reps)]
    bad = [r for r in recs if not r["ok"]]
    if bad:
        sys.exit(f"reference run failed on {n}; see {R.logdir}/{bad[0].get('log')}")
    ref = dict(B=recs[0]["B"], TF=recs[0]["TF"], TFSm=recs[0]["TFSm"],
               NB=recs[0]["NB"], M=recs[0]["M"],
               t_total=statistics.median(r["t_total"] for r in recs),
               wall=statistics.median(r["wall"] for r in recs),
               spread=(max(r["t_total"] for r in recs) - min(r["t_total"] for r in recs)) /
               statistics.mean(r["t_total"] for r in recs),
               reps=reps, cpu=recs[0].get("cpu"))
    meta["refs"][n] = ref
    return ref


# --------------------------------------------------------------------------
# prune lists
# --------------------------------------------------------------------------
def prune_path(outdir, d):
    return os.path.join(outdir, f"prune_c{d}.txt")


def load_prune(path):
    levels, points = [], set()
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            kind, rest = line.split(None, 1)
            if kind == "level":
                k, v = rest.replace(" ", "").split("=")
                levels.append((k, float(v)))
            elif kind == "point":
                points.add(point_key(parse_point(rest)))
            else:
                raise ValueError(f"{path}: bad line {line!r}")
    return levels, points


def is_pruned(pt, prune):
    if prune is None:
        return False
    levels, points = prune
    if point_key(pt) in points:
        return True
    return any(abs(float(pt[k]) - v) < 1e-9 for k, v in levels)


# --------------------------------------------------------------------------
# collecting results from the cache
# --------------------------------------------------------------------------
def collect(outdir, meta, inputs=None):
    """one entry per (point, number) run: ratio = t_total / reference time.
    Timed-out runs are 'censored' (slower than timeout_mult x reference)."""
    cache = v1.load_cache(outdir)
    out = []
    for idx, n in enumerate(meta["numbers"], 1):
        if (inputs and idx not in inputs) or n not in meta["refs"]:
            continue
        ref = meta["refs"][n]
        for pt in grid_points(meta["levels"]):
            k = v1.run_key(n, abs_params(pt, ref, meta), 0, meta["threads"], meta["extra"])
            rec = cache.get(k)
            if rec is None:
                continue
            e = dict(point=pt, pkey=point_key(pt), idx=idx, n=n, rec=rec,
                     censored=bool(rec.get("timeout")),
                     failed=(not rec["ok"]) and not rec.get("timeout"))
            e["ratio"] = rec["t_total"] / ref["t_total"] if rec["ok"] else None
            out.append(e)
    return out


def point_scores(data, timeout_mult):
    """per point: geometric mean ratio over the numbers it ran on (censored
    runs count as timeout_mult, failed runs exclude the point)"""
    by = {}
    for e in data:
        by.setdefault(e["pkey"], []).append(e)
    scores = {}
    for pk, es in by.items():
        if any(e["failed"] for e in es):
            scores[pk] = dict(point=es[0]["point"], score=float("inf"), count=len(es),
                              censored=False, failed=True, ratios=[])
            continue
        rs = [timeout_mult if e["censored"] else e["ratio"] for e in es]
        scores[pk] = dict(point=es[0]["point"], score=v1.gmean(rs), count=len(es),
                          censored=any(e["censored"] for e in es), failed=False,
                          ratios=rs)
    return scores


# --------------------------------------------------------------------------
# grid
# --------------------------------------------------------------------------
def parse_inputs(s, count):
    if s in (None, "all"):
        return list(range(1, count + 1))
    return [int(x) for x in s.split(",") if x]


def levels_from_args(a, d):
    lv = default_levels(d)
    for k in GRID_PARAMS:
        v = getattr(a, "lv_" + k)
        if v:
            lv[k] = [float(x) if k == "B" else int(x) for x in v.split(",")]
    return lv


def cmd_grid(a):
    check_binary(a)
    R = make_runner(a)
    for d in a.digits:
        meta = init_meta(a, d)
        levels = levels_from_args(a, d)
        for k in GRID_PARAMS:       # the meta keeps the union of levels ever run
            meta["levels"][k] = sorted(set(meta["levels"][k]) | set(levels[k]))
        save_meta(a.outdir, meta)
        prune = None
        if a.prune:
            pp = prune_path(a.outdir, d) if a.prune == "auto" else a.prune.format(d=d)
            prune = load_prune(pp)
            log(f"c{d}: pruning with {pp}: {len(prune[0])} levels, {len(prune[1])} points")

        for idx in parse_inputs(a.inputs, len(meta["numbers"])):
            n = meta["numbers"][idx - 1]
            ref = reference(R, meta, n, a.ref_reps, a.max_time)
            save_meta(a.outdir, meta)
            timeout = max(a.min_timeout, a.timeout_mult * ref["wall"])
            pts = [pt for pt in grid_points(levels) if not is_pruned(pt, prune)]
            random.Random(f"{d}:{idx}:{a.seed}").shuffle(pts)
            log(f"\n===== c{d} input {idx}: {n} =====")
            log(f"reference (-noopt, {v1.fmt_params(meta['pipeline'])}): "
                f"B={ref['B']} TF={ref['TF']} TFSm(F)={ref['TFSm']} NB={ref['NB']} "
                f"M={ref['M']}  t={ref['t_total']:.2f}s (spread {100 * ref['spread']:.1f}%)")
            log(f"{len(pts)} points, timeout {timeout:.1f}s")
            t0 = time.time()
            best = None
            for i, pt in enumerate(pts, 1):
                rec = R.run(n, abs_params(pt, ref, meta), 0, timeout)
                if rec.get("timeout"):
                    res = f">{a.timeout_mult:.2f} (timeout)"
                elif not rec["ok"]:
                    res = f"FAILED (log {rec.get('log')})"
                else:
                    r = rec["t_total"] / ref["t_total"]
                    if best is None or r < best[0]:
                        best = (r, point_key(pt))
                    res = f"{r:.3f}"
                eta = (time.time() - t0) / i * (len(pts) - i)
                bs = f"best {best[0]:.3f} ({best[1]})" if best else ""
                log(f"[c{d} #{idx} {i}/{len(pts)}] {point_key(pt):40s} {res:18s} "
                    f"{bs}  ETA {eta / 60:.0f} min")
        save_meta(a.outdir, meta)
    log(f"\n{R.nrun} new yafu runs; cache {R.cache_path}")


# --------------------------------------------------------------------------
# analyze
# --------------------------------------------------------------------------
def cmd_analyze(a):
    for d in a.digits:
        meta = load_meta(a.outdir, d)
        inputs = parse_inputs(a.inputs, len(meta["numbers"])) if a.inputs else None
        data = collect(a.outdir, meta, inputs)
        if not data:
            log(f"c{d}: no grid runs yet")
            continue
        used = sorted({e["idx"] for e in data})
        scores = point_scores(data, a.timeout_mult)
        good = [s for s in scores.values() if not s["failed"]]
        best = min(good, key=lambda s: s["score"])
        nfail = sum(s["failed"] for s in scores.values())
        ncens = sum(e["censored"] for e in data)
        print(f"\n===== c{d}: {len(data)} runs on input(s) {used}, {len(scores)} points "
              f"({ncens} timed out, {nfail} points failed) =====")
        for idx in used:
            n = meta["numbers"][idx - 1]
            r = meta["refs"][n]
            print(f"  input {idx}: reference B={r['B']} TF={r['TF']} F={r['TFSm']} "
                  f"NB={r['NB']} M={r['M']} t={r['t_total']:.2f}s")
        print("ratio = time / reference (yafu defaults, -noopt, fixed pipeline); "
              "geometric mean over the inputs a point ran on")

        print("\nper level: best / median / worst point ratio")
        for k in GRID_PARAMS:
            for v in meta["levels"][k]:
                ss = [s["score"] for s in good if abs(float(s["point"][k]) - v) < 1e-9]
                if not ss:
                    continue
                print(f"  {k:4s} {fmt_level(k, v):>5s}:  {min(ss):.3f} / "
                      f"{statistics.median(ss):.3f} / {max(ss):.3f}   ({len(ss)} points)")

        for k1, k2 in (("B", "M"), ("TF", "TFSm"), ("B", "NB"), ("B", "TF")):
            v1s, v2s = meta["levels"][k1], meta["levels"][k2]
            print(f"\nbest ratio by {k1} (rows) x {k2} (columns)")
            print("  " + " " * 7 + "".join(f"{fmt_level(k2, v):>8s}" for v in v2s))
            for x in v1s:
                cells = []
                for y in v2s:
                    ss = [s["score"] for s in good
                          if abs(float(s["point"][k1]) - x) < 1e-9 and
                          abs(float(s["point"][k2]) - y) < 1e-9]
                    cells.append(f"{min(ss):8.3f}" if ss else f"{'-':>8s}")
                print(f"  {fmt_level(k1, x):>7s}" + "".join(cells))

        print(f"\nbest {a.top} points (x = inputs measured on)")
        for s in sorted(good, key=lambda s: s["score"])[:a.top]:
            print(f"  {s['score']:.3f} x{s['count']}  {point_key(s['point']):40s} "
                  f"({' '.join(f'{r:.3f}' for r in s['ratios'])})")

        # prune suggestions
        b = best["score"]
        lines = [f"# siqs_tune_v2 prune list for c{d}, from input(s) {used}; best "
                 f"{b:.3f} ({point_key(best['point'])})",
                 f"# 'level' lines drop a whole level whose best point is > "
                 f"{100 * a.level_tol:.0f}% slower than the best; 'point' lines "
                 f"drop points > {100 * a.point_tol:.0f}% slower (or timed out / "
                 f"failed).  Delete lines to keep them."]
        pruned_levels = []
        for k in GRID_PARAMS:
            for v in meta["levels"][k]:
                ss = [s["score"] for s in good if abs(float(s["point"][k]) - v) < 1e-9]
                if ss and min(ss) > b * (1 + a.level_tol):
                    lines.append(f"level {k}={fmt_level(k, v)}   # best {min(ss):.3f}")
                    pruned_levels.append((k, v))
        npts = 0
        for s in sorted(scores.values(), key=lambda s: s["score"]):
            pt = s["point"]
            if any(abs(float(pt[k]) - v) < 1e-9 for k, v in pruned_levels):
                continue
            if s["failed"] or s["censored"] or s["score"] > b * (1 + a.point_tol):
                why = "failed" if s["failed"] else ("timed out" if s["censored"]
                                                    else f"{s['score']:.3f}")
                lines.append(f"point {point_key(pt)}   # {why}")
                npts += 1
        total = len(grid_points(meta["levels"]))
        keep = [pt for pt in grid_points(meta["levels"])
                if not is_pruned(pt, (pruned_levels, {l.split()[1] for l in lines
                                                      if l.startswith("point")}))]
        print(f"\nsuggested prune: {len(pruned_levels)} levels, {npts} more points -> "
              f"{len(keep)} of {total} points left per further input")
        if not a.no_write:
            pp = prune_path(a.outdir, d)
            with open(pp, "w") as f:
                f.write("\n".join(lines) + "\n")
            print(f"wrote {pp} (review it, then: grid --inputs 2,3 --prune)")


# --------------------------------------------------------------------------
# fit: quadratic response surface
# --------------------------------------------------------------------------
def coords(pt, centers):
    """scaled coordinates: b = log2(B/0.7), m = log2(M/180), t, s, nb centered"""
    return [math.log2(float(pt["B"]) / 0.7), math.log2(pt["M"] / 180.0),
            (pt["TF"] - centers["TF"]) / 4.0, (pt["TFSm"] - centers["TFSm"]) / 2.0,
            (pt["NB"] - centers["NB"])]


def make_terms(nlevels):
    """the surface's terms, given how many distinct levels each coordinate has
    in the fit data: none for a single level, linear for two, linear and
    square for three or more; cross terms between every pair that varies.
    A term is a tuple of coordinate indices; () is the intercept."""
    act = [i for i in range(FIT_TERMS) if nlevels[i] >= 2]
    terms = [()] + [(i,) for i in act] + [(i, i) for i in act if nlevels[i] >= 3]
    terms += [(i, j) for n, i in enumerate(act) for j in act[n + 1:]]
    return terms


def features(x, terms):
    return [math.prod(x[i] for i in t) for t in terms]


def term_name(t):
    names = ["b", "m", "t", "s", "nb"]
    return "1" if not t else "*".join(names[i] for i in t)


def solve(A, y, ridge=1e-6):
    """least squares via normal equations (small systems; pure python)"""
    p = len(A[0])
    N = [[sum(r[i] * r[j] for r in A) + (ridge if i == j and i > 0 else 0.0)
          for j in range(p)] for i in range(p)]
    v = [sum(r[i] * yy for r, yy in zip(A, y)) for i in range(p)]
    for c in range(p):                      # gaussian elimination, partial pivot
        piv = max(range(c, p), key=lambda r: abs(N[r][c]))
        N[c], N[piv] = N[piv], N[c]
        v[c], v[piv] = v[piv], v[c]
        if abs(N[c][c]) < 1e-12:
            continue
        for r in range(c + 1, p):
            f = N[r][c] / N[c][c]
            if f:
                for k in range(c, p):
                    N[r][k] -= f * N[c][k]
                v[r] -= f * v[c]
    beta = [0.0] * p
    for c in range(p - 1, -1, -1):
        if abs(N[c][c]) < 1e-12:
            continue
        beta[c] = (v[c] - sum(N[c][k] * beta[k] for k in range(c + 1, p))) / N[c][c]
    return beta


def cmd_fit(a):
    for d in a.digits:
        meta = load_meta(a.outdir, d)
        data = [e for e in collect(a.outdir, meta) if e["ratio"] is not None]
        if not data:
            log(f"c{d}: no completed grid runs")
            continue
        lv = meta["levels"]
        centers = {k: statistics.median(lv[k]) for k in ("TF", "TFSm", "NB")}
        idxs = sorted({e["idx"] for e in data})
        # fit near the optimum: each number's runs within fit_within of its best
        bestn = {i: min(e["ratio"] for e in data if e["idx"] == i) for i in idxs}
        sel = [e for e in data if e["ratio"] <= bestn[e["idx"]] * (1 + a.fit_within)]
        nlevels = [len({float(e["point"][k]) for e in sel}) for k in GRID_PARAMS]
        terms = make_terms(nlevels)
        nfeat = len(terms)
        A, y = [], []
        for e in sel:
            row = features(coords(e["point"], centers), terms)
            row += [1.0 if e["idx"] == i else 0.0 for i in idxs[1:]]   # per-number intercepts
            A.append(row)
            y.append(math.log(e["ratio"]))
        need = nfeat + len(idxs) - 1 + 5
        if len(A) < need:
            log(f"c{d}: only {len(A)} runs within {100 * a.fit_within:.0f}% of the best, "
                f"{need} needed for {nfeat} terms ({' '.join(map(term_name, terms))}); "
                f"widen --fit-within or run more inputs")
            continue
        beta = solve(A, y)
        pred = [sum(b * x for b, x in zip(beta, r)) for r in A]
        ybar = statistics.mean(y)
        ss_res = sum((yy - pp) ** 2 for yy, pp in zip(y, pred))
        ss_tot = sum((yy - ybar) ** 2 for yy in y) or 1e-12
        rmse = math.sqrt(ss_res / max(1, len(y) - len(beta)))
        # intercept averaged over the numbers, so predictions are "typical input"
        off = statistics.mean([0.0] + beta[nfeat:])

        def predict(pt):
            return math.exp(sum(b * x for b, x in zip(beta[:nfeat],
                                                      features(coords(pt, centers), terms))) + off)

        # optimum over the surveyed box (no extrapolation)
        span = {k: (min(float(e["point"][k]) for e in sel), max(float(e["point"][k]) for e in sel))
                for k in GRID_PARAMS}
        nb_steps = max(1, int(round((span["B"][1] - span["B"][0]) / 0.02)))
        Bs = [span["B"][0] + i * (span["B"][1] - span["B"][0]) / nb_steps
              for i in range(nb_steps + 1)]
        Ms = [int(round(span["M"][0] * (span["M"][1] / span["M"][0]) ** (i / 12.0)))
              for i in range(13)] if span["M"][1] > span["M"][0] else [int(span["M"][0])]
        rngi = lambda k: list(range(int(span[k][0]), int(span[k][1]) + 1))  # noqa: E731
        best_pt, best_val = None, float("inf")
        for B, M, TF, TS, NB in itertools.product(Bs, sorted(set(Ms)), rngi("TF"),
                                                  rngi("TFSm"), rngi("NB")):
            pt = dict(B=round(B, 3), M=M, TF=TF, TFSm=TS, NB=NB)
            v = predict(pt)
            if v < best_val:
                best_pt, best_val = pt, v

        sens = {}
        steps = dict(B=("x", 1.15), M=("x", 1.4), TF=("+", 2), TFSm=("+", 2), NB=("+", 1))
        for k, (kind, s) in steps.items():
            vals = []
            for sign in (-1, 1):
                pt = dict(best_pt)
                if kind == "x":
                    pt[k] = pt[k] * (s if sign > 0 else 1 / s)
                else:
                    pt[k] = pt[k] + sign * s
                inside = span[k][0] - 1e-9 <= float(pt[k]) <= span[k][1] + 1e-9
                vals.append(predict(pt) / best_val if inside else None)
            sens[k] = vals

        scores = point_scores(collect(a.outdir, meta), a.timeout_mult)
        full = [s for s in scores.values() if s["count"] == len(idxs) and not s["failed"]]
        bm = min(full, key=lambda s: s["score"]) if full else None
        refs = [meta["refs"][meta["numbers"][i - 1]] for i in idxs]
        refmed = {k: statistics.median(r[k] for r in refs) for k in ("B", "TF", "TFSm", "NB", "M")}
        absopt = dict(B=v1.quantize("B", best_pt["B"] * refmed["B"]), M=best_pt["M"],
                      TF=int(refmed["TF"] + best_pt["TF"]),
                      TFSm=int(refmed["TFSm"] + best_pt["TFSm"]), NB=best_pt["NB"])

        res = dict(digits=d, inputs=idxs, n_obs=len(sel), n_runs=len(data),
                   fit_within=a.fit_within, r2=1 - ss_res / ss_tot, rmse=rmse,
                   optimum=best_pt, pred_ratio=best_val, sensitivity=sens,
                   optimum_abs=absopt, ref_median=refmed,
                   best_measured=dict(point=bm["point"], score=bm["score"],
                                      count=bm["count"]) if bm else None,
                   span=span, terms=[term_name(t) for t in terms], beta=beta)
        with open(os.path.join(a.outdir, f"fit_c{d}.json"), "w") as f:
            json.dump(res, f, indent=1)

        print(f"\n===== c{d}: fit on {len(sel)} of {len(data)} runs (each input's runs "
              f"within {100 * a.fit_within:.0f}% of its best), inputs {idxs} =====")
        print(f"{nfeat} terms: {' '.join(map(term_name, terms))}  (b=log2 B/0.7, "
              f"m=log2 M/180, t=TF, s=TFSm, nb=NB; one level = left out, two = linear)")
        print(f"R^2 = {res['r2']:.3f}, residual sd = {100 * rmse:.1f}% (run noise ~1-2%)")
        print(f"fitted optimum: {point_key(best_pt)}  predicted ratio {best_val:.3f}")
        print(f"  = {v1.fmt_params(absopt)} for a typical c{d} (reference medians "
              f"B={refmed['B']:.0f} TF={refmed['TF']:.0f} F={refmed['TFSm']:.0f})")
        if bm:
            print(f"best measured on all {len(idxs)} input(s): {point_key(bm['point'])} "
                  f"{bm['score']:.3f}")
        print("sensitivity (predicted ratio vs optimum, one step down / up):")
        for k, (kind, s) in steps.items():
            lo, hi = sens[k]
            step = f"/{s}, x{s}" if kind == "x" else f"-{s}, +{s}"
            print(f"  {k:4s} ({step:>12s}): "
                  f"{(f'{lo:.3f}' if lo else '  -  '):>6s}  {(f'{hi:.3f}' if hi else '  -  '):>6s}")


# --------------------------------------------------------------------------
# trends across sizes
# --------------------------------------------------------------------------
def linfit(xs, ys):
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return 0.0, my
    sl = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return sl, my - sl * mx


def cmd_trends(a):
    fits = []
    for p in os.listdir(a.outdir):
        if p.startswith("fit_c") and p.endswith(".json"):
            with open(os.path.join(a.outdir, p)) as f:
                fits.append(json.load(f))
    if not fits:
        sys.exit(f"no fit_c*.json in {a.outdir}: run 'fit' first")
    fits.sort(key=lambda r: r["digits"])
    print(f"{'size':>5s} {'inputs':>7s} {'B/def':>6s} {'M':>5s} {'TF':>4s} {'TFSm':>5s} "
          f"{'NB':>3s} {'pred':>6s} {'R2':>5s}   absolute optimum            best measured")
    for r in fits:
        o = r["optimum"]
        bm = r.get("best_measured")
        print(f"c{r['digits']:<4d} {len(r['inputs']):>7d} {o['B']:6.2f} {o['M']:5d} "
              f"{o['TF']:+4d} {o['TFSm']:+5d} {o['NB']:3d} {r['pred_ratio']:6.3f} "
              f"{r['r2']:5.2f}   {v1.fmt_params(r['optimum_abs']):28s} "
              + (f"{point_key(bm['point'])} {bm['score']:.3f}" if bm else ""))
    if len(fits) >= 2:
        xs = [r["digits"] for r in fits]
        print("\nlinear trends vs digits (TF offset from yafu's default TF, TFSm "
              "from the formula F):")
        for k, get in (("B/def", lambda r: r["optimum"]["B"]),
                       ("log2 M", lambda r: math.log2(r["optimum"]["M"])),
                       ("TF", lambda r: r["optimum"]["TF"]),
                       ("TFSm", lambda r: r["optimum"]["TFSm"]),
                       ("NB", lambda r: r["optimum"]["NB"])):
            sl, ic = linfit(xs, [get(r) for r in fits])
            preds = "  ".join(f"c{p}: {(2 ** (sl * p + ic) if k == 'log2 M' else sl * p + ic):.2f}"
                              for p in a.predict)
            print(f"  {k:7s} slope {sl:+.4f}/digit   predicted {preds}")
        print("(use these to center narrower grids for larger sizes: grid --B ... --TF ...)")


# --------------------------------------------------------------------------
# refine (pass 2) -- writes siqs_tune.py-compatible summaries
# --------------------------------------------------------------------------
def cmd_refine(a):
    check_binary(a)
    R = make_runner(a)
    for d in a.digits:
        meta = load_meta(a.outdir, d)
        if (meta["threads"], meta["extra"]) != (a.threads, a.extra):
            sys.exit(f"c{d}: refine with the grid's --threads {meta['threads']} "
                     f"--extra '{meta['extra']}'")
        train = [n for n in meta["numbers"] if n in meta["refs"]]
        refs = {n: meta["refs"][n] for n in train}
        noise = max(0.015, statistics.mean(r["spread"] for r in refs.values()))
        timeout = max(a.min_timeout, a.timeout_mult * max(r["wall"] for r in refs.values()))
        rng = random.Random(f"refine:{d}:{a.seed}")

        # center: the fit's optimum, the best measured grid point, or explicit
        refmed = {k: statistics.median(r[k] for r in refs.values())
                  for k in ("B", "TF", "TFSm", "NB")}
        if a.center in ("fit", "best"):
            fp = os.path.join(a.outdir, f"fit_c{d}.json")
            if a.center == "fit" and os.path.exists(fp):
                with open(fp) as f:
                    pt = json.load(f)["optimum"]
            else:
                sc = point_scores(collect(a.outdir, meta), a.timeout_mult)
                full = [s for s in sc.values() if s["count"] == len(train) and not s["failed"]]
                if not full:
                    sys.exit(f"c{d}: no grid point measured on all inputs; run fit or grid first")
                pt = min(full, key=lambda s: s["score"])["point"]
                if a.center == "fit":
                    log(f"c{d}: no fit_c{d}.json, centering on the best measured point")
            center = dict(B=v1.quantize("B", pt["B"] * refmed["B"]), M=int(pt["M"]),
                          TF=int(refmed["TF"] + pt["TF"]),
                          TFSm=int(refmed["TFSm"] + pt["TFSm"]), NB=int(pt["NB"]))
        else:
            center = {k: v1.quantize(k, v if k in v1.STR_CHOICES else float(v))
                      for k, v in v1.parse_kv(a.center).items()}
        center.setdefault("MFBD", meta["mfbd"])
        for k, v in meta["pipeline"].items():
            center.setdefault(k, v)

        names = ["B", "M", "NB", "TF", "TFSm"]
        fixed = {}
        if d >= 70:
            names.append("MFBD")
            if a.free_pipeline:
                names += ["Mid", "Large", "Kernel"]
        for k in ("MFBD", "Mid", "Large", "Kernel"):
            if k not in names:
                fixed[k] = center[k]
        box = {k: (0, 0) for k in v1.STR_CHOICES}
        box.update(B=(v1.quantize("B", center["B"] * 0.8), v1.quantize("B", center["B"] * 1.25)),
                   M=(max(20, int(center["M"] * 0.7)), int(center["M"] * 1.4) + 1),
                   NB=(max(1, center["NB"] - 1), center["NB"] + 1),
                   TF=(center["TF"] - 3, center["TF"] + 3),
                   TFSm=(max(1, center["TFSm"] - 3), center["TFSm"] + 3),
                   MFBD=(round(center["MFBD"] - 0.1, 2), round(center["MFBD"] + 0.1, 2)))
        log(f"\n===== c{d} refine around {v1.fmt_params(center)} on {len(train)} "
            f"grid inputs, {a.trials} trials =====")
        log("box: " + "  ".join(f"{k}=[{box[k][0]},{box[k][1]}]" for k in names
                                if k not in v1.STR_CHOICES) +
            (f"  fixed: {v1.fmt_params(fixed)}" if fixed else ""))

        cands = {}

        def evaluate(c, rep):
            ratios = []
            for n in train:
                rec = R.run(n, c.params, rep, timeout)
                if not rec["ok"]:
                    c.failed = True
                    return
                ratios.append(rec["t_total"] / refs[n]["t_total"])
            c.reps.append(v1.gmean(ratios))

        def incumbent():
            pool = [c for c in cands.values() if len(c.reps) >= 2] or list(cands.values())
            return min(pool, key=lambda c: c.score)

        t0 = time.time()

        def consider(params, origin, label):
            params = dict(params)
            k = v1.pkey(params)
            if k in cands:
                return
            c = v1.Candidate(params, origin)
            cands[k] = c
            evaluate(c, 0)
            best = min((x.score for x in cands.values() if len(x.reps) >= 2), default=1.0)
            r = 1
            while not c.failed and r < a.max_reps and c.score < best * (1 + 2 * noise):
                evaluate(c, r)
                r += 1
            inc = incumbent()
            log(f"[c{d} {label} {origin:6s}] {v1.fmt_params(params):60s} "
                f"{'FAIL' if c.failed else f'{c.score:.3f}'} x{len(c.reps)}  "
                f"best {inc.score:.3f}  {(time.time() - t0) / 60:.1f} min")

        consider(center, "center", "0")
        for t in range(a.trials):
            prog = t / max(1, a.trials)
            scale = 0.3 * (1 - prog) + 0.08
            consider(v1.perturb(incumbent().params, names, box, fixed, scale, rng),
                     "local", f"{t + 1}/{a.trials}")

        # holdout check, siqs_tune.py style
        ranked = sorted((c for c in cands.values() if not c.failed), key=lambda c: c.score)
        by_key = {v1.pkey(c.params): c for c in ranked}
        picked = v1.select_diverse([(c.params, c.score) for c in ranked], a.final_top,
                                   a.diverse_tol, [k for k in names if k not in fixed], box)
        top = [by_key[v1.pkey(p)] for p, _ in picked]
        hold = v1.test_numbers(d, a.holdout, a.seed, "holdout2")
        log(f"\nverifying {len(top)} sets on {len(hold)} holdout numbers, "
            f"{a.final_reps} interleaved rounds")
        contenders = [("asis", {}, None), ("builtin", v1.BUILTIN_BUCKET, None)] + \
            [("tuned", c.params, c) for c in top]
        times = {i: {n: [] for n in hold} for i in range(len(contenders))}
        hrecs = {i: {n: [] for n in hold} for i in range(len(contenders))}
        cpu = None
        for rep in range(a.final_reps):
            for n in hold:
                for i, (kind, p, _) in enumerate(contenders):
                    rec = R.run(n, p, 100 + rep,
                                timeout=timeout if kind == "tuned" else a.max_time)
                    times[i][n].append(rec["t_total"] if rec["ok"] else float("inf"))
                    hrecs[i][n].append(rec)
                    cpu = cpu or rec.get("cpu")
        results = []
        for i, (kind, p, c) in enumerate(contenders):
            per_n = [statistics.median(times[i][n]) / statistics.median(times[0][n]) for n in hold]
            sp = v1.gmean(per_n) if all(math.isfinite(x) for x in per_n) else float("inf")
            se = statistics.stdev(per_n) / math.sqrt(len(per_n)) if len(per_n) > 1 else 0.0
            results.append(dict(params=p, kind=kind, ratio=sp, se=se,
                                train=c.score if c else (1.0 if kind == "asis" else None),
                                mean_time=statistics.mean(statistics.median(times[i][n])
                                                          for n in hold),
                                stats=v1.run_stats(hrecs[i])))
        r0 = refs[train[0]]
        summary = dict(digits=d, bits=statistics.mean(int(n).bit_length() for n in hold),
                       cpu=cpu, names=names,
                       defaults=dict(B=r0["B"], M=r0["M"], NB=r0["NB"], TF=r0["TF"],
                                     TFSm=r0["TFSm"], **meta["pipeline"]),
                       base_times={n: refs[n]["t_total"] for n in train}, noise=noise,
                       box={k: box[k] for k in names}, fixed=fixed, train=train,
                       holdout=hold, results=results, objective="total",
                       threads=a.threads, extra=a.extra, source="siqs_tune_v2 refine",
                       candidates=[dict(params=c.params, origin=c.origin, reps=c.reps,
                                        failed=c.failed) for c in ranked])
        with open(os.path.join(a.outdir, f"summary_c{d}.json"), "w") as f:
            json.dump(summary, f, indent=1)
        v1.print_summary(summary)
    log(f"\n{R.nrun} new yafu runs; cache {R.cache_path}")


def cmd_export(a):
    v1.cmd_export(argparse.Namespace(outdir=a.outdir, top=a.top, max_gap=a.max_gap,
                                     numbers=a.numbers, out=a.out))


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, runs=False):
        p.add_argument("--outdir", default="siqs_tune_v2_out")
        p.add_argument("--digits", type=int, nargs="+", default=[70])
        p.add_argument("--timeout-mult", type=float, default=1.5,
                       help="runs slower than this x the reference are stopped "
                            "(and count as this in analyze); default 1.5")
        if runs:
            p.add_argument("--yafu", required=True)
            p.add_argument("--threads", type=int, default=1)
            p.add_argument("--extra", default="", help="extra yafu args for every run")
            p.add_argument("--ini", help="yafu.ini to copy into each work dir")
            p.add_argument("--keep-logs", action="store_true")
            p.add_argument("--seed", default="1")
            p.add_argument("--min-timeout", type=float, default=10.0,
                           help="never stop a run sooner than this (seconds)")
            p.add_argument("--max-time", type=float, default=3600,
                           help="timeout for reference / default runs")
            p.add_argument("--allow-binary-change", action="store_true")

    g = sub.add_parser("grid", help="pass 1: nested sweep")
    common(g, runs=True)
    g.add_argument("--inputs", default="1", help="which test numbers: 1, 2,3 or all")
    g.add_argument("--numbers", type=int, default=3, help="test numbers per size")
    g.add_argument("--ref-reps", type=int, default=2, help="reference repeats per number")
    g.add_argument("--pipeline", default=DEFAULT_PIPELINE,
                   help=f"fixed Mid/Large/Kernel (default {DEFAULT_PIPELINE})")
    g.add_argument("--mfbd", type=float, default=1.72)
    g.add_argument("--prune", nargs="?", const="auto",
                   help="skip points in a prune list (default <outdir>/prune_cNN.txt; "
                        "a path may use {d} for the size)")
    for k in GRID_PARAMS:
        g.add_argument(f"--{k}", dest="lv_" + k,
                       help=f"comma list of {k} levels (overrides the size default)")
    g.set_defaults(func=cmd_grid)

    an = sub.add_parser("analyze", help="tables + suggested prune list")
    common(an)
    an.add_argument("--inputs", help="restrict to these inputs (default all run)")
    an.add_argument("--top", type=int, default=15)
    an.add_argument("--level-tol", type=float, default=0.08,
                    help="prune a level whose best point is this much slower than the best")
    an.add_argument("--point-tol", type=float, default=0.15,
                    help="prune points this much slower than the best")
    an.add_argument("--no-write", action="store_true", help="don't write prune_cNN.txt")
    an.set_defaults(func=cmd_analyze)

    fi = sub.add_parser("fit", help="quadratic response surface per size")
    common(fi)
    fi.add_argument("--fit-within", type=float, default=0.25,
                    help="fit each input's runs within this fraction of its best")
    fi.set_defaults(func=cmd_fit)

    tr = sub.add_parser("trends", help="fitted optima across sizes")
    tr.add_argument("--outdir", default="siqs_tune_v2_out")
    tr.add_argument("--predict", type=int, nargs="+", default=[85, 90, 95])
    tr.set_defaults(func=cmd_trends)

    rf = sub.add_parser("refine", help="pass 2: local search + holdout check")
    common(rf, runs=True)
    rf.add_argument("--center", default="fit",
                    help="fit (default), best (best measured grid point), or "
                         "explicit absolute params B=..,M=..,NB=..,TF=..,TFSm=..")
    rf.add_argument("--trials", type=int, default=40)
    rf.add_argument("--max-reps", type=int, default=3)
    rf.add_argument("--free-pipeline", action="store_true",
                    help="also search Mid/Large/Kernel (>= c70)")
    rf.add_argument("--holdout", type=int, default=3)
    rf.add_argument("--final-top", type=int, default=5)
    rf.add_argument("--final-reps", type=int, default=2)
    rf.add_argument("--diverse-tol", type=float, default=0.05)
    rf.set_defaults(func=cmd_refine)

    ex = sub.add_parser("export", help="siqs_candidates.txt (siqs_tune.py's export)")
    ex.add_argument("--outdir", nargs="+", default=["siqs_tune_v2_out"])
    ex.add_argument("--top", type=int, default=5)
    ex.add_argument("--max-gap", type=float, default=0.10)
    ex.add_argument("--numbers", type=int, default=1)
    ex.add_argument("--out")
    ex.set_defaults(func=cmd_export)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
