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
  SPV   -siqsSPV   small prime variation bound: primes below it are not
                   sieved (one of 256, 512, 1024; default 256).  Only has an
                   effect in builds with TRY_SMALLP_OPT; yafu ignores the
                   flag otherwise, so drop it with --params on such builds.
  Mid   -siqsMid   sieve for primes 32k..interval: sparse or bucket
  Large -siqsLarge sieve for primes > interval: direct or bucket
  Kernel -siqsKernel direct-sieve scatter kernel: cmem, creg, sink, bitloop
                   (Mid/Large/Kernel are categorical; their defaults are read
                   from yafu's "strategy:" line.  Builds without the
                   large-prime strategies don't print it, and the three are
                   then dropped.)

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

  # write the competing parameter sets for yafu's tune (per digit size: the
  # best, the most different near-best sets, and yafu's defaults), with the
  # partial relations each needed to finish -> siqs_candidates.txt / .json.
  # Copy siqs_candidates.txt next to yafu and run siqstune() there: it picks
  # the fastest set per size on that machine and writes siqs_params.txt.
  python3 siqs_tune.py export --outdir siqs_tune_out

Note: the direct sieve needs a DLP-sized input (>= 70 digits); below that
yafu falls back to large-prime buckets, so Large=direct is a no-op there.

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
    "SPV":      (re.compile(r"SPV small_limit = (\d+)"), int),   # -v -v only
    "Mid":      (re.compile(r"strategy: mid = (\w+)"), str),
    "Large":    (re.compile(r"strategy: mid = \w+, large = (\w+)"), str),
    "Kernel":   (re.compile(r"strategy: .*?, kernel = (\w+)"), str),
    # results
    "t_total":  (re.compile(r"SIQS elapsed time = " + FLT), float),
    "t_sieve":  (re.compile(r"^QS elapsed time = " + FLT, re.M), float),
    "t_lanczos": (re.compile(r"Lanczos elapsed time = " + FLT), float),
    "polys":    (re.compile(r"sieving required (\d+) total polynomials"), int),
    "a_polys":  (re.compile(r"total polynomials \((\d+) 'A'"), int),
    "rels_sec": (re.compile(r"\(" + FLT + r" rels/sec\)"), float),
    "needed":   (re.compile(r"(\d+) relations needed"), int),
    # "N rels found: F full + C from P partial" (screen progress line) or
    # "N relations found: ..." (final summary); the last match wins
    "partials": (re.compile(r"(?:rels|relations) found:\s*\d+ full \+ \d+ from (\d+) partial"), int),
    "fulls":    (re.compile(r"(?:rels|relations) found:\s*(\d+) full \+"), int),
    "cpu":      (re.compile(r"^Detected (?!L1 )(.+?)\s*$", re.M), str),
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
         "TF": "-siqsTF", "TFSm": "-siqsTFSm", "MFBD": "-siqsMFBD",
         "SPV": "-siqsSPV", "Mid": "-siqsMid", "Large": "-siqsLarge",
         "Kernel": "-siqsKernel"}
ORDER = ["B", "M", "NB", "TF", "TFSm", "MFBD", "SPV", "Mid", "Large", "Kernel"]
LOG_PARAMS = {"B", "M"}
# categorical params: the box (lo, hi) selects the allowed choices
CHOICES = {"SPV": [256, 512, 1024]}
SPV_DEFAULT = 256           # yafu's small_limit; not printed at -v
# string-valued categorical params (the large-prime strategy); no box
STR_CHOICES = {"Mid": ["sparse", "bucket"],
               "Large": ["direct", "bucket"],
               "Kernel": ["cmem", "creg", "sink", "bitloop"]}


def choices(name, box):
    if name in STR_CHOICES:
        return STR_CHOICES[name]
    lo, hi = box[name]
    return [c for c in CHOICES[name] if lo <= c <= hi] or [lo]


def is_choice(name):
    return name in CHOICES or name in STR_CHOICES


def default_box(d):
    """Search ranges around yafu's defaults d (dict of parsed effective params)."""
    box = {
        "B":    (max(256, int(d["B"] * 0.35)), int(d["B"] * 1.3)),
        "M":    (20, 600),
        "NB":   (1, max(8, d["NB"] + 2)),
        "TF":   (max(10, d["TF"] - 14), d["TF"] + 3),
        # a raised SPV bound leaves more small primes unsieved, so the
        # matching TFSm correction can sit a few bits above the default's
        "TFSm": (max(6, d["TFSm"] - 8), d["TFSm"] + 6),
        "MFBD": (1.45, 1.95),
        "SPV":  (min(CHOICES["SPV"]), max(CHOICES["SPV"])),
    }
    for k in STR_CHOICES:
        box[k] = (0, 0)             # unused; choices come from STR_CHOICES
    return box


def quantize(name, v):
    if name in STR_CHOICES:
        return str(v)
    if name == "B":
        return max(16, int(round(v / 16.0)) * 16)
    if name == "MFBD":
        return round(float(v), 2)
    return int(round(v))


def to_unit(name, v, box):
    if name in STR_CHOICES:
        cs = choices(name, box)
        i = cs.index(v) if v in cs else 0
        return i / (len(cs) - 1) if len(cs) > 1 else 0.5
    if name in CHOICES:
        cs = choices(name, box)
        i = min(range(len(cs)), key=lambda j: abs(cs[j] - v))
        return i / (len(cs) - 1) if len(cs) > 1 else 0.5
    lo, hi = box[name]
    if name in LOG_PARAMS:
        return (math.log(v) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return (v - lo) / (hi - lo) if hi > lo else 0.5


def from_unit(name, u, box):
    lo, hi = box[name]
    u = min(1.0, max(0.0, u))
    if is_choice(name):
        cs = choices(name, box)
        return cs[min(len(cs) - 1, int(u * len(cs)))]
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
        if is_choice(k):             # no notion of step size: jump to another
            others = [c for c in choices(k, box) if c != center[k]]
            if others:
                p[k] = rng.choice(others)
            continue
        for _ in range(20):          # make sure the quantized value changes
            v = from_unit(k, to_unit(k, center[k], box) + rng.gauss(0, scale), box)
            if v != center[k]:
                break
        p[k] = v
    return p


def pkey(params):
    return ",".join(f"{k}={params[k]}" for k in ORDER if k in params)


def run_key(n, params, rep, threads, extra):
    return f"{n}|{pkey(params)}|{rep}|t{threads}|{extra}"


def load_cache(outdir):
    cache = {}
    path = os.path.join(outdir, "runs.jsonl")
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue
                cache[rec["key"]] = rec
    return cache


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
        self.cache = load_cache(a.outdir)
        self.nrun = 0

    def key(self, n, params, rep):
        return run_key(n, params, rep, self.a.threads, self.a.extra)

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
        if "supported options are:" in text:
            # yafu prints its usage and exits on an option it doesn't know;
            # don't cache this as a failed candidate
            raise SystemExit(f"{self.yafu} rejected the command line (printed usage):\n"
                             f"  {' '.join(cmd)}\n"
                             "is the binary older than one of these options? "
                             "(drop it with --params)")
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


STAT_FIELDS = ("partials", "fulls", "needed", "t_sieve", "t_post", "t_total",
               "rels_sec")


def run_stats(recs_by_n):
    """Per-parameter-set statistics for the yafu tune table.

    recs_by_n: {n: [records]}.  Each number's repeats are reduced to their
    median first, then summarized across numbers: the median, plus the
    coefficient of variation of partials (how well 'partials needed'
    transfers between inputs of the same size)."""
    per = {k: [] for k in STAT_FIELDS}
    for n, recs in recs_by_n.items():
        ok = [r for r in recs if r["ok"]]
        if not ok:
            return None
        for k in STAT_FIELDS:
            if k == "t_post":
                vs = [r["t_total"] - r["t_sieve"] for r in ok
                      if r.get("t_sieve") is not None]
            else:
                vs = [r[k] for r in ok if r.get(k) is not None]
            if vs:
                per[k].append(statistics.median(vs))
    out = {k: (statistics.median(vs) if vs else None) for k, vs in per.items()}
    ps = per["partials"]
    out["partials_cv"] = (statistics.stdev(ps) / statistics.mean(ps)
                          if len(ps) > 1 and statistics.mean(ps) > 0 else None)
    return out


def pdist(p, q, names, box):
    """Distance between parameter sets in the search's unit space (log for B
    and M); a differing string choice (Mid/Large/Kernel) counts as 1."""
    s = 0.0
    for k in names:
        if k not in p or k not in q:
            continue
        if k in STR_CHOICES:
            s += float(p[k] != q[k])
        else:
            s += (to_unit(k, p[k], box) - to_unit(k, q[k], box)) ** 2
    return math.sqrt(s / max(1, len(names)))


def select_diverse(ranked, k, tol, names, box, alpha=0.25):
    """Pick k parameter sets: the best, then greedily the near-best sets
    (score within tol of the best) farthest from those already picked, with
    a mild penalty for being slower.  If fewer than k sets are near-best,
    fill with the next best by score.

    ranked: list of (params, score) sorted by score."""
    if not ranked:
        return []
    best = ranked[0][1]
    pool = [x for x in ranked if x[1] <= best * (1 + tol)]
    chosen = [pool.pop(0)]
    while pool and len(chosen) < k:
        def value(x):
            dmin = min(pdist(x[0], c[0], names, box) for c in chosen)
            return dmin - alpha * (x[1] / best - 1) / tol
        x = max(pool, key=value)
        pool.remove(x)
        chosen.append(x)
    for x in ranked:
        if len(chosen) >= k:
            break
        if all(x is not c for c in chosen):
            chosen.append(x)
    return chosen


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
                eff.setdefault("SPV", SPV_DEFAULT)
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
        if k in STR_CHOICES:
            sys.exit(f"--range doesn't apply to {k}; use --fix {k}=<value>")
        lo, hi = rng_s.split(":")
        box[k] = (quantize(k, float(lo)), quantize(k, float(hi)))
    fixed = {}
    for k, v in parse_kv(",".join(a.fix or [])).items():
        fixed[k] = eff[k] if v == "default" else \
            (quantize(k, v) if k in STR_CHOICES else quantize(k, float(v)))
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
    # verify the best plus the near-best sets that differ most from it: those
    # are the ones other CPUs may rank differently
    by_key = {pkey(c.params): c for c in ranked}
    picked = select_diverse([(c.params, c.score) for c in ranked
                             if c.origin != "default"],
                            a.final_top, a.diverse_tol,
                            [k for k in names if k not in fixed], box)
    top = [by_key[pkey(p)] for p, _ in picked]
    log(f"\nverifying {len(top)} sets (best + most different within "
        f"{100 * a.diverse_tol:.0f}% of it) on {len(hold)} holdout numbers, "
        f"{a.final_reps} interleaved rounds")
    contenders = [("defaults", {})] + [(fmt_params(c.params), c.params) for c in top]
    times = {i: {n: [] for n in hold} for i in range(len(contenders))}
    hrecs = {i: {n: [] for n in hold} for i in range(len(contenders))}
    cpu = None
    for rep in range(a.final_reps):
        for n in hold:
            for i, (_, p) in enumerate(contenders):
                # holdout baseline is yafu's own default (no flags)
                rec = R.run(n, p, 100 + rep, timeout=timeout if p else a.max_time)
                t = metric(rec, a.objective)
                times[i][n].append(t if t is not None else float("inf"))
                hrecs[i][n].append(rec)
                cpu = cpu or rec.get("cpu")

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
                                statistics.median(times[i][n]) for n in hold),
                            stats=run_stats(hrecs[i])))
    bits = statistics.mean(int(n).bit_length() for n in hold)
    summary = dict(digits=digits, bits=bits, cpu=cpu, names=names,
                   defaults=eff, base_times=base_t, noise=noise,
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
    print(f"{'holdout t/t_def':>16s} {'+-se':>6s} {'train':>6s} {'mean s':>7s} "
          f"{'partials':>9s} {'cv':>5s}  params")
    for r in s["results"]:
        label = fmt_params(r["params"]) if r["params"] else "(yafu defaults)"
        st = r.get("stats") or {}
        pt = f"{st['partials']:9.0f}" if st.get("partials") else f"{'-':>9s}"
        cv = (f"{st['partials_cv']:5.2f}" if st.get("partials_cv") is not None
              else f"{'-':>5s}")
        print(f"{r['ratio']:16.3f} {r['se']:6.3f} {r['train']:6.3f} "
              f"{r['mean_time']:7.2f} {pt} {cv}  {label}")
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


def describe_diff(p, ref):
    """Short description of how parameter set p differs from ref."""
    out = []
    for k in ORDER:
        if k not in p or k not in ref or p[k] == ref[k]:
            continue
        if k in LOG_PARAMS:
            out.append(f"{k} x{p[k] / ref[k]:.2f}")
        elif k in STR_CHOICES or k in CHOICES:
            out.append(f"{k}={p[k]}")
        elif k == "MFBD":
            out.append(f"{k} {p[k] - ref[k]:+.2f}")
        else:
            out.append(f"{k} {p[k] - ref[k]:+d}")
    return ", ".join(out) or "same"


EXPORT_COLS = ["digits", "bits", "rank", "ratio"] + ORDER + \
    ["partials", "partials_cv", "fulls", "needed", "t_sieve", "t_post", "n"]


def cmd_export(a):
    """Write the competing parameter sets for yafu's tune: per digit size, the
    holdout-verified sets (best first, then the most different near-best
    ones) plus yafu's defaults as rank 0.

    Each set gets one row per exported holdout number n, with the relation
    counts and times measured on that n: yafu's tune re-runs the sets on the
    same n, where the partials needed to finish reproduce to ~0.5%."""
    # one summary per digit size; with several outdirs, a later one wins
    summaries = {}
    for d in a.outdir:
        for p in os.listdir(d):
            if p.startswith("summary_c") and p.endswith(".json"):
                digits = int(re.findall(r"\d+", p)[0])
                if digits in summaries:
                    log(f"c{digits}: using {d}, not {summaries[digits][0]}")
                summaries[digits] = (d, p)
    if not summaries:
        sys.exit(f"no summary_c*.json in {', '.join(a.outdir)}")
    caches = {}
    rows, cpus = [], set()
    for digits in sorted(summaries):
        outdir, p = summaries[digits]
        with open(os.path.join(outdir, p)) as f:
            s = json.load(f)
        if s.get("cpu"):
            cpus.add(s["cpu"])
        if outdir not in caches:
            caches[outdir] = load_cache(outdir)
        cache = caches[outdir]
        suffix = f"|t{s['threads']}|{s['extra']}"

        def holdout_recs(params, n):
            return [rec for rec in cache.values()
                    if rec["n"] == n and rec["rep"] >= 100 and
                    rec["key"].endswith(suffix) and
                    pkey(rec["params"]) == pkey(params)]

        for r in s["results"]:
            if r.get("stats") is None:
                # summary written before stats were recorded: rebuild them
                # from the holdout runs in the cache
                recs = {n: holdout_recs(r["params"], n) for n in s["holdout"]}
                r["stats"] = run_stats(recs) if all(recs.values()) else None
        res = [r for r in s["results"] if r["params"] and math.isfinite(r["ratio"])]
        res.sort(key=lambda r: r["ratio"])
        if not res:
            continue
        best = res[0]
        keep = [best] + [r for r in res[1:] if r["ratio"] <= best["ratio"] * (1 + a.max_gap)]
        keep = keep[:a.top]
        default = next(r for r in s["results"] if not r["params"])
        entries = [(0, {}, default)] + \
                  [(i + 1, r["params"], r) for i, r in enumerate(keep)]
        numbers = s["holdout"][:a.numbers] if a.numbers > 0 else s["holdout"]
        for n in numbers:
            for rank, params, r in entries:
                recs = [x for x in holdout_recs(params, n) if x["ok"]]
                if not recs:
                    log(f"c{digits} rank {rank}: no holdout run on {n[:12]}..., skipped")
                    continue
                st = run_stats({n: recs}) or {}
                if rank == 0:
                    # run without flags; list the values yafu chose for this n
                    params = {k: recs[0][k] for k in ORDER if recs[0].get(k) is not None}
                row = dict(digits=s["digits"], bits=int(n).bit_length(), rank=rank,
                           ratio=round(r["ratio"], 4), n=n)
                row.update({k: params.get(k) for k in ORDER})
                row.update({k: st.get(k) for k in EXPORT_COLS if k in st})
                row["vs_best"] = describe_diff(params, best["params"]) if rank > 1 else ""
                rows.append(row)

    base = a.out or os.path.join(a.outdir[0], "siqs_candidates")
    with open(base + ".json", "w") as f:
        json.dump(dict(cpu=sorted(cpus), rows=rows), f, indent=1)

    def cell(v):
        if v is None:
            return "-"
        if isinstance(v, float):
            return f"{v:.3f}" if abs(v) < 100 else f"{v:.0f}"
        return str(v)

    with open(base + ".txt", "w") as f:
        f.write(f"# siqs_tune export; cpu: {'; '.join(sorted(cpus)) or '?'}\n")
        f.write("# rank 0 = yafu defaults (run without flags; params are the "
                "values yafu chose for this n), rank 1 = best tuned set on "
                "this cpu\n")
        f.write("# one row per set per test number n; partials, fulls, "
                "needed, t_sieve, t_post = measured on that n (median over "
                "repeats); ratio = holdout time / defaults time on this cpu, "
                "over all holdout numbers\n")
        f.write("# " + " ".join(EXPORT_COLS) + "\n")
        for row in rows:
            f.write(" ".join(cell(row.get(k)) for k in EXPORT_COLS) + "\n")

    # human-readable view
    last = None
    for row in rows:
        if (row["digits"], row["n"]) != last:
            last = (row["digits"], row["n"])
            print(f"\n===== c{row['digits']} ({row['bits']} bits): {row['n']} =====")
            print(f"{'rank':>4s} {'ratio':>6s} {'partials':>9s} {'cv':>5s} "
                  f"{'t_post':>6s}  params / difference from best")
        pv = row.get("partials")
        cv = row.get("partials_cv")
        tp = row.get("t_post")
        print(f"{row['rank']:4d} {row['ratio']:6.3f} "
              f"{(f'{pv:9.0f}' if pv else '-'):>9s} "
              f"{(f'{cv:5.2f}' if cv is not None else '-'):>5s} "
              f"{(f'{tp:6.2f}' if tp is not None else '-'):>6s}  "
              f"{fmt_params({k: row[k] for k in ORDER if row.get(k) is not None})}"
              + (f"\n{'':36s}({row['vs_best']})" if row["vs_best"] else
                 ("  (yafu defaults)" if row["rank"] == 0 else "")))
    print(f"\nwrote {base}.txt and {base}.json")


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
    t.add_argument("--final-top", type=int, default=5,
                   help="sets verified on holdout: the best + most different near-best")
    t.add_argument("--diverse-tol", type=float, default=0.05,
                   help="near-best = train score within this fraction of the best")
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

    e = sub.add_parser("export", help="write the competing parameter sets "
                       "(table for yafu's tune)")
    e.add_argument("--outdir", nargs="+", default=["siqs_tune_out"],
                   help="tune output dir(s); several dirs are merged")
    e.add_argument("--top", type=int, default=5,
                   help="max sets per digit size, best included (default 5)")
    e.add_argument("--max-gap", type=float, default=0.10,
                   help="drop sets more than this fraction slower than the "
                        "best on the holdout (default 0.10)")
    e.add_argument("--numbers", type=int, default=1,
                   help="holdout numbers to export per digit size, each with "
                        "its own rows (default 1; 0 = all).  yafu's tune "
                        "benchmarks every set on each")
    e.add_argument("--out", help="output path prefix (default "
                   "<first outdir>/siqs_candidates -> .txt and .json)")
    e.set_defaults(func=cmd_export)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
