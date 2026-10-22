#!/usr/bin/env python3
"""
siqs_tune.py -- search for better YAFU SIQS parameters at a given input size.

For each requested digit size the tuner:
  1. builds a deterministic set of balanced semiprimes (the "train" set) and a
     second, disjoint set (the "holdout" set) used only for final verification;
  2. runs yafu with its own defaults on the train set ('baseline'), and parses
     the parameters yafu actually chose (factor base size, LP multiplier,
     blocks, TF bound, small-prime-variation cutoff, DLP MFB);
  3. builds a search box around those defaults and runs a Monte-Carlo /
     local-perturbation search.  Each candidate is scored by the geometric
     mean over the train numbers of  t_candidate / t_baseline,  so < 1.0 is
     faster.  Candidates that look competitive get re-run (racing) so a lucky
     single run can't win;
  4. re-runs the top few candidates and the defaults, interleaved, on the
     holdout set and reports speedups with a noise estimate.

Every yafu run is cached in <outdir>/runs.jsonl keyed by (number, params,
repeat), so an interrupted tune resumes where it stopped and a re-run with
more --trials only does the new work.

Tunable parameters (yafu flag):
  B     -siqsB     factor base size, primes      (log-uniform, multiple of 16)
  M     -siqsM     large prime multiplier        (log-uniform)
  NB    -siqsNB    number of 32k sieve blocks
  TF    -siqsTF    trial-factoring cutoff, bits
  TFSm  -siqsTFSm  small prime variation cutoff (also disables yafu's
                   adaptive TFSm tuning, as any -siqsTFSm override does)
  MFBD  -siqsMFBD  DLP residue exponent (only meaningful when DLP is used,
                   i.e. >= 70 digits; dropped automatically otherwise)

Examples
  # tune at 70 digits, 1 thread, 60 candidates
  python3 siqs_tune.py tune --yafu ./yafu --digits 70 --trials 60

  # also seed the search with a point known to be good
  python3 siqs_tune.py tune --yafu ./yafu --digits 73 \\
      --seed-point "B=11008,M=300,NB=3,TF=75,TFSm=15,MFBD=1.67"

  # restrict / widen ranges, hold a parameter fixed
  python3 siqs_tune.py tune --yafu ./yafu --digits 70 \\
      --range NB=1:6 --range M=50:600 --fix TFSm=default

  # summarize what's in an output directory
  python3 siqs_tune.py report --outdir siqs_tune_out

Note: builds with TRY_DIRECT_SIEVE exit below the DLP cutoff (70 digits),
so tune at >= 70 digits with those builds.

yafu treats a non-tty stdin as a batchfile, so on Linux/WSL each run gets a
pseudo-terminal; on Windows stdin is inherited from the console.
"""

import argparse
import concurrent.futures as cf
import hashlib
import json
import math
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time

try:
    import pty
except ImportError:          # Windows
    pty = None


# --------------------------------------------------------------------------
# test number generation (pure python, deterministic)
# --------------------------------------------------------------------------
SMALL_PRIMES = [p for p in range(3, 2000)
                if all(p % q for q in range(2, int(p ** 0.5) + 1))]


def is_prp(n, rounds=24, rng=None):
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n % p == 0:
            return n == p
    rng = rng or random.Random(n)
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def rand_prime(rng, lo, hi):
    while True:
        p = rng.randrange(lo, hi) | 1
        if is_prp(p, rng=rng):
            return p


def semiprime(digits, rng):
    """Balanced semiprime with exactly `digits` decimal digits."""
    lo, hi = 10 ** (digits - 1), 10 ** digits
    while True:
        p = rand_prime(rng, math.isqrt(lo) + 1, math.isqrt(hi))
        qlo, qhi = max(lo // p + 1, p // 2), min(hi // p, 2 * p)
        if qlo >= qhi:
            continue
        q = rand_prime(rng, qlo, qhi)
        n = p * q
        if lo <= n < hi and p != q:
            return n


def test_numbers(digits, count, seed, kind):
    # string seed -> stable across python versions/runs
    rng = random.Random(f"{kind}:{digits}:{seed}")
    return [str(semiprime(digits, rng)) for _ in range(count)]


# --------------------------------------------------------------------------
# yafu output parsing
# --------------------------------------------------------------------------
FLT = r"([0-9]*\.?[0-9]+)"
PARSE = {
    # effective parameters
    "B":        (re.compile(r"factor base:\s*(\d+) primes"), int),
    "pmax":     (re.compile(r"factor base:.*max prime = (\d+)"), int),
    "M":        (re.compile(r"single large prime cutoff:\s*\d+\s*\((\d+)\s*\*\s*pmax"), int),
    "MFBD":     (re.compile(r"DLP MFB = " + FLT), float),
    "NB":       (re.compile(r"sieve interval:\s*(\d+) blocks"), int),
    "TFSm":     (re.compile(r"using SPV correction of (\d+) bits"), int),
    "TF":       (re.compile(r"trial factoring cutoff at (\d+) bits"), int),
    # results
    "t_total":  (re.compile(r"SIQS elapsed time = " + FLT), float),
    "t_sieve":  (re.compile(r"^QS elapsed time = " + FLT, re.M), float),
    "t_lanczos": (re.compile(r"Lanczos elapsed time = " + FLT), float),
    "polys":    (re.compile(r"sieving required (\d+) total polynomials"), int),
    "a_polys":  (re.compile(r"total polynomials \((\d+) 'A'"), int),
    "rels_sec": (re.compile(r"\(" + FLT + r" rels/sec\)"), float),
    "needed":   (re.compile(r"(\d+) relations needed"), int),
}
RE_FACTORS = re.compile(r"\*\*\*factors found\*\*\*")
RE_FACTOR = re.compile(r"^[PC]\d+\s*=\s*(\d+)", re.M)


def parse_output(text, n):
    r = {}
    for k, (rx, conv) in PARSE.items():
        ms = list(rx.finditer(text))
        r[k] = conv(ms[-1].group(1)) if ms else None
    ok = False
    if RE_FACTORS.search(text):
        fs = [int(f) for f in RE_FACTOR.findall(text.split("***factors found***")[-1])]
        nn = int(n)
        ok = len(fs) >= 2 and math.prod(fs) == nn and all(1 < f < nn for f in fs)
    r["ok"] = ok and r["t_total"] is not None
    return r


# --------------------------------------------------------------------------
# parameter space
# --------------------------------------------------------------------------
FLAGS = {"B": "-siqsB", "M": "-siqsM", "NB": "-siqsNB",
         "TF": "-siqsTF", "TFSm": "-siqsTFSm", "MFBD": "-siqsMFBD"}
ORDER = ["B", "M", "NB", "TF", "TFSm", "MFBD"]
LOG_PARAMS = {"B", "M"}


def default_box(d):
    """Search ranges around yafu's defaults d (dict of parsed effective params)."""
    box = {
        "B":    (max(256, int(d["B"] * 0.35)), int(d["B"] * 1.3)),
        "M":    (20, 600),
        "NB":   (1, max(8, d["NB"] + 2)),
        "TF":   (max(10, d["TF"] - 14), d["TF"] + 3),
        "TFSm": (max(6, d["TFSm"] - 8), d["TFSm"] + 3),
        "MFBD": (1.45, 1.95),
    }
    return box


def quantize(name, v):
    if name == "B":
        return max(16, int(round(v / 16.0)) * 16)
    if name == "MFBD":
        return round(float(v), 2)
    return int(round(v))


def to_unit(name, v, box):
    lo, hi = box[name]
    if name in LOG_PARAMS:
        return (math.log(v) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return (v - lo) / (hi - lo) if hi > lo else 0.5


def from_unit(name, u, box):
    lo, hi = box[name]
    u = min(1.0, max(0.0, u))
    if name in LOG_PARAMS:
        v = math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo)))
    else:
        v = lo + u * (hi - lo)
    return quantize(name, v)


def sample_uniform(names, box, fixed, rng):
    p = {k: from_unit(k, rng.random(), box) for k in names}
    p.update(fixed)
    return p


def perturb(center, names, box, fixed, scale, rng):
    """Gaussian step in unit space on a random subset (>=1) of the free params."""
    free = [k for k in names if k not in fixed]
    if not free:
        return dict(center)
    move = [k for k in free if rng.random() < 0.5] or [rng.choice(free)]
    p = dict(center)
    for k in move:
        for _ in range(20):          # make sure the quantized value changes
            v = from_unit(k, to_unit(k, center[k], box) + rng.gauss(0, scale), box)
            if v != center[k]:
                break
        p[k] = v
    return p


def pkey(params):
    return ",".join(f"{k}={params[k]}" for k in ORDER if k in params)


def parse_kv(s):
    out = {}
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        k, v = tok.split("=", 1)
        k = k.strip()
        if k not in FLAGS:
            sys.exit(f"unknown parameter '{k}' (known: {', '.join(ORDER)})")
        out[k] = v.strip()
    return out


# --------------------------------------------------------------------------
# running yafu, with a persistent cache
# --------------------------------------------------------------------------
class Runner:
    def __init__(self, a):
        self.a = a
        self.yafu = os.path.abspath(a.yafu)
        self.ini = os.path.abspath(a.ini) if a.ini else None
        os.makedirs(a.outdir, exist_ok=True)
        self.cache_path = os.path.join(a.outdir, "runs.jsonl")
        self.logdir = os.path.join(a.outdir, "logs")
        if a.keep_logs:
            os.makedirs(self.logdir, exist_ok=True)
        self.lock = threading.Lock()
        self.cache = {}
        self.nrun = 0
        if os.path.exists(self.cache_path):
            with open(self.cache_path) as f:
                for line in f:
                    try:
                        rec = json.loads(line)
                    except ValueError:
                        continue
                    self.cache[rec["key"]] = rec

    def key(self, n, params, rep):
        return f"{n}|{pkey(params)}|{rep}|t{self.a.threads}|{self.a.extra}"

    def cmd(self, n, params):
        c = [self.yafu, f"siqs({n})", "-v", "-threads", str(self.a.threads),
             "-inmem", "1000"]
        for k in ORDER:
            if k in params:
                c += [FLAGS[k], str(params[k])]
        if self.a.extra:
            c += self.a.extra.split()
        return c

    def run(self, n, params, rep, timeout):
        k = self.key(n, params, rep)
        with self.lock:
            if k in self.cache:
                return self.cache[k]
        tag = hashlib.sha1(k.encode()).hexdigest()[:12]
        wd = os.path.join(self.a.outdir, "work", tag)
        shutil.rmtree(wd, ignore_errors=True)
        os.makedirs(wd)
        if self.ini:
            shutil.copy(self.ini, os.path.join(wd, "yafu.ini"))
        cmd = self.cmd(n, params)
        t0 = time.time()
        timed_out = False
        master = slave = None
        try:
            if pty is not None:
                master, slave = pty.openpty()
                stdin = slave
            else:
                stdin = None
            p = subprocess.run(cmd, cwd=wd, stdin=stdin, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, timeout=timeout)
            text = p.stdout.decode(errors="replace")
        except subprocess.TimeoutExpired as e:
            timed_out = True
            text = (e.stdout or b"").decode(errors="replace") + "\n### TIMEOUT\n"
        finally:
            if slave is not None:
                os.close(slave)
                os.close(master)
        wall = time.time() - t0
        shutil.rmtree(wd, ignore_errors=True)
        rec = parse_output(text, n)
        rec.update(key=k, n=n, digits=len(n), params=params, rep=rep,
                   wall=round(wall, 3), timeout=timed_out, when=time.time())
        if self.a.keep_logs or not rec["ok"]:
            os.makedirs(self.logdir, exist_ok=True)
            with open(os.path.join(self.logdir, tag + ".log"), "w") as f:
                f.write("### " + " ".join(cmd) + "\n" + text)
            rec["log"] = tag + ".log"
        with self.lock:
            self.cache[k] = rec
            self.nrun += 1
            with open(self.cache_path, "a") as f:
                f.write(json.dumps(rec) + "\n")
        return rec

    def run_many(self, jobs, timeout):
        """jobs: list of (n, params, rep).  Returns records in the same order."""
        if self.a.jobs <= 1 or len(jobs) == 1:
            return [self.run(n, p, r, timeout) for n, p, r in jobs]
        with cf.ThreadPoolExecutor(self.a.jobs) as ex:
            return list(ex.map(lambda j: self.run(j[0], j[1], j[2], timeout), jobs))


def metric(rec, objective):
    if not rec["ok"]:
        return None
    if objective == "sieve":
        return rec["t_sieve"] or rec["t_total"]
    return rec["t_total"]


def gmean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


# --------------------------------------------------------------------------
# tuning
# --------------------------------------------------------------------------
class Candidate:
    def __init__(self, params, origin):
        self.params = params
        self.origin = origin
        self.reps = []          # list of per-repeat score (gmean ratio)
        self.failed = False

    @property
    def score(self):
        if self.failed or not self.reps:
            return float("inf")
        return statistics.mean(self.reps)


def log(*args):
    print(*args, file=sys.stderr, flush=True)


def fmt_params(p):
    return " ".join(f"{k}={p[k]}" for k in ORDER if k in p)


def yafu_flags(p):
    return " ".join(f"{FLAGS[k]} {p[k]}" for k in ORDER if k in p)


def tune_digits(R, a, digits):
    rng = random.Random(f"search:{digits}:{a.seed}")
    train = test_numbers(digits, a.numbers, a.seed, "train")
    hold = test_numbers(digits, a.holdout, a.seed, "holdout")
    log(f"\n===== c{digits}: {len(train)} train numbers, {len(hold)} holdout =====")

    # ---- baseline -------------------------------------------------------
    base = {}                   # n -> list of times
    eff = None
    for rep in range(a.base_reps):
        recs = R.run_many([(n, {}, rep) for n in train], timeout=a.max_time)
        for n, rec in zip(train, recs):
            t = metric(rec, a.objective)
            if t is None:
                sys.exit(f"baseline run failed on {n}; see {R.logdir}/{rec.get('log')}")
            base.setdefault(n, []).append(t)
            if eff is None:
                eff = {k: rec[k] for k in ORDER if rec.get(k) is not None}
    base_t = {n: statistics.median(ts) for n, ts in base.items()}
    noise = [statistics.stdev(ts) / statistics.mean(ts) for ts in base.values() if len(ts) > 1]
    noise = statistics.mean(noise) if noise else 0.03
    log(f"defaults chosen by yafu: {fmt_params(eff)}")
    log("baseline times: " + ", ".join(f"{base_t[n]:.2f}" for n in train) +
        f"  (run-to-run noise ~{100 * noise:.1f}%)")
    timeout = max(30.0, a.timeout_mult * max(base_t.values()))

    # ---- search space -----------------------------------------------------
    names = [k for k in ORDER if k in eff]          # MFBD absent below DLP cutoff
    if a.params:
        names = [k for k in names if k in a.params.split(",")]
    box = default_box(eff)
    for spec in a.range or []:
        k, rng_s = spec.split("=", 1)
        lo, hi = rng_s.split(":")
        box[k] = (quantize(k, float(lo)), quantize(k, float(hi)))
    fixed = {}
    for k, v in parse_kv(",".join(a.fix or [])).items():
        fixed[k] = eff[k] if v == "default" else quantize(k, float(v))
        if k not in names:
            names.append(k)
    log("search box: " + "  ".join(f"{k}=[{box[k][0]},{box[k][1]}]"
                                   for k in names if k not in fixed) +
        ("  fixed: " + fmt_params(fixed) if fixed else ""))

    cands = {}

    def evaluate(c, rep):
        recs = R.run_many([(n, c.params, rep) for n in train], timeout=timeout)
        ratios = []
        for n, rec in zip(train, recs):
            t = metric(rec, a.objective)
            if t is None:
                c.failed = True
                return
            ratios.append(t / base_t[n])
        c.reps.append(gmean(ratios))

    t_start = time.time()

    def consider(params, origin, label):
        params = {k: params[k] for k in names}
        k = pkey(params)
        if k in cands:
            return None
        c = Candidate(params, origin)
        cands[k] = c
        evaluate(c, 0)
        # racing: confirm anything close to the incumbent before trusting it
        best = min((x.score for x in cands.values() if len(x.reps) >= 2),
                   default=1.0)
        r = 1
        while (not c.failed and r < a.max_reps and
               c.score < best * (1 + 2 * noise)):
            evaluate(c, r)
            r += 1
        inc = incumbent()
        log(f"[c{digits} {label} {origin:7s}] {fmt_params(c.params):48s} "
            f"{'FAIL' if c.failed else f'{c.score:.3f}'} x{len(c.reps)}   "
            f"best {inc.score:.3f} ({fmt_params(inc.params)})   "
            f"{(time.time() - t_start) / 60:.1f} min")
        return c

    def incumbent():
        pool = [c for c in cands.values() if len(c.reps) >= 2] or list(cands.values())
        return min(pool, key=lambda c: c.score)

    # the defaults themselves, expressed as explicit flags, are candidate 0;
    # this also measures the cost of disabling adaptive TFSm.
    consider(dict(eff), "default", "0")
    for sp in a.seed_point or []:
        p = dict(eff)
        p.update({k: quantize(k, float(v)) for k, v in parse_kv(sp).items()})
        p.update(fixed)
        consider(p, "seed", "0")

    n_explore = int(a.trials * a.explore)
    for t in range(a.trials):
        if t < n_explore:
            params, origin = sample_uniform(names, box, fixed, rng), "random"
        else:
            prog = (t - n_explore) / max(1, a.trials - n_explore)
            scale = 0.25 * (1 - prog) + 0.04
            params = perturb(incumbent().params, names, box, fixed, scale, rng)
            origin = "local"
        consider(params, origin, f"{t + 1}/{a.trials}")

    # ---- holdout verification --------------------------------------------
    ranked = sorted((c for c in cands.values() if not c.failed),
                    key=lambda c: c.score)
    top = []
    for c in ranked:
        if len(top) >= a.final_top:
            break
        if c.origin == "default":
            continue
        top.append(c)
    log(f"\nverifying top {len(top)} on {len(hold)} holdout numbers, "
        f"{a.final_reps} interleaved rounds")
    contenders = [("defaults", {})] + [(fmt_params(c.params), c.params) for c in top]
    times = {i: {n: [] for n in hold} for i in range(len(contenders))}
    for rep in range(a.final_reps):
        for n in hold:
            for i, (_, p) in enumerate(contenders):
                # holdout baseline is yafu's own default (no flags)
                rec = R.run(n, p, 100 + rep, timeout=timeout if p else a.max_time)
                t = metric(rec, a.objective)
                times[i][n].append(t if t is not None else float("inf"))

    results = []
    for i, (label, p) in enumerate(contenders):
        per_n = []
        for n in hold:
            bt = statistics.median(times[0][n])
            per_n.append(statistics.median(times[i][n]) / bt)
        sp = gmean(per_n) if all(math.isfinite(x) for x in per_n) else float("inf")
        spread = statistics.stdev(per_n) / math.sqrt(len(per_n)) if len(per_n) > 1 else 0.0
        results.append(dict(params=p, ratio=sp, se=spread,
                            train=top[i - 1].score if i else 1.0,
                            mean_time=statistics.mean(
                                statistics.median(times[i][n]) for n in hold)))
    summary = dict(digits=digits, defaults=eff, base_times=base_t, noise=noise,
                   box={k: box[k] for k in names}, fixed=fixed, train=train,
                   holdout=hold, results=results, objective=a.objective,
                   threads=a.threads, extra=a.extra,
                   candidates=[dict(params=c.params, origin=c.origin, reps=c.reps,
                                    failed=c.failed) for c in ranked])
    with open(os.path.join(a.outdir, f"summary_c{digits}.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print_summary(summary)
    return summary


def print_summary(s):
    print(f"\n===== c{s['digits']}  (objective: {s['objective']} time, "
          f"threads={s['threads']}) =====")
    print(f"yafu defaults: {fmt_params(s['defaults'])}")
    print(f"{'holdout t/t_def':>16s} {'+-se':>6s} {'train':>6s} {'mean s':>7s}  params")
    for r in s["results"]:
        label = fmt_params(r["params"]) if r["params"] else "(yafu defaults)"
        print(f"{r['ratio']:16.3f} {r['se']:6.3f} {r['train']:6.3f} "
              f"{r['mean_time']:7.2f}  {label}")
    best = min(s["results"], key=lambda r: r["ratio"])
    if best["params"]:
        print(f"best: {100 * (1 - best['ratio']):.1f}% faster than defaults on holdout"
              f" -> {yafu_flags(best['params'])}")
    else:
        print("no candidate beat the defaults on the holdout set")


def cmd_tune(a):
    R = Runner(a)
    for d in a.digits:
        tune_digits(R, a, d)
    log(f"\n{R.nrun} new yafu runs; cache: {R.cache_path}")


def cmd_report(a):
    paths = sorted(p for p in os.listdir(a.outdir)
                   if p.startswith("summary_c") and p.endswith(".json"))
    if not paths:
        sys.exit(f"no summary_c*.json in {a.outdir}")
    for p in sorted(paths, key=lambda s: int(re.findall(r"\d+", s)[0])):
        with open(os.path.join(a.outdir, p)) as f:
            s = json.load(f)
        print_summary(s)
        if a.top:
            print(f"top {a.top} train candidates:")
            for c in s["candidates"][:a.top]:
                sc = statistics.mean(c["reps"]) if c["reps"] else float("nan")
                print(f"  {sc:.3f} x{len(c['reps'])} {c['origin']:7s} {fmt_params(c['params'])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("tune", help="run the parameter search")
    t.add_argument("--yafu", required=True, help="path to yafu binary")
    t.add_argument("--digits", type=int, nargs="+", default=[70],
                   help="input sizes to tune (default 70)")
    t.add_argument("--outdir", default="siqs_tune_out")
    t.add_argument("--threads", type=int, default=1, help="yafu -threads")
    t.add_argument("--jobs", type=int, default=1,
                   help="concurrent yafu runs (shares L2/L3/turbo: biases timing)")
    t.add_argument("--numbers", type=int, default=3, help="train numbers per size")
    t.add_argument("--holdout", type=int, default=3, help="holdout numbers per size")
    t.add_argument("--seed", default="1", help="seed for numbers and search")
    t.add_argument("--trials", type=int, default=60, help="candidates to try")
    t.add_argument("--explore", type=float, default=0.4,
                   help="fraction of trials sampled uniformly before local search")
    t.add_argument("--base-reps", type=int, default=3,
                   help="baseline repeats per train number")
    t.add_argument("--max-reps", type=int, default=3,
                   help="max repeats for a competitive candidate (racing)")
    t.add_argument("--final-top", type=int, default=4,
                   help="candidates verified on holdout")
    t.add_argument("--final-reps", type=int, default=2,
                   help="interleaved rounds in holdout verification")
    t.add_argument("--objective", choices=["total", "sieve"], default="total",
                   help="'total' = SIQS elapsed incl. filtering+LA (default); "
                        "'sieve' = sieving only")
    t.add_argument("--params", help="comma list of params to tune (default all)")
    t.add_argument("--range", action="append", metavar="P=LO:HI",
                   help="override a search range, e.g. NB=1:6 (repeatable)")
    t.add_argument("--fix", action="append", metavar="P=V",
                   help="hold a param at V or 'default', e.g. TFSm=default")
    t.add_argument("--seed-point", action="append", metavar="P=V,...",
                   help="extra starting candidate, unspecified params = default")
    t.add_argument("--extra", default="", help="extra yafu args for every run")
    t.add_argument("--ini", help="yafu.ini to copy into each run's work dir")
    t.add_argument("--timeout-mult", type=float, default=3.0,
                   help="kill candidates slower than this x slowest baseline")
    t.add_argument("--max-time", type=float, default=3600,
                   help="timeout for default-parameter runs, seconds")
    t.add_argument("--keep-logs", action="store_true",
                   help="keep every yafu log (failed runs are always kept)")
    t.set_defaults(func=cmd_tune)

    r = sub.add_parser("report", help="print summaries from an outdir")
    r.add_argument("--outdir", default="siqs_tune_out")
    r.add_argument("--top", type=int, default=10, help="also list N best train candidates")
    r.set_defaults(func=cmd_report)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
