#!/usr/bin/env python3
"""
Empirical test of YAFU's ecm() factor-finding ability.

For each of N trials:
  1. build n = p*q with p ~ B/3 bits and q ~ 2B/3 bits (random primes)
  2. run:  yafu "ecm(<n>, C)" -B1ecm <B1>
  3. parse the portion of factor.log written by that run
After all trials, print summary statistics (factors found, curves run, ...).

Example:
  python3 ecm_test.py -B 150 -N 50 -C 200 --B1 11000 --yafu ./yafu
"""
import argparse
import os
import random
import re
import shlex
import statistics
import subprocess
import sys
import time

# ----------------------------------------------------------------- primes
_SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]


def is_prime(n, rng, rounds=40):
    if n < 2:
        return False
    for sp in _SMALL:
        if n == sp:
            return True
        if n % sp == 0:
            return False
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


def random_prime(bits, rng):
    """Random prime with exactly `bits` bits."""
    if bits < 2:
        raise ValueError("prime size must be >= 2 bits")
    while True:
        c = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(c, rng):
            return c


def make_input(B, rng):
    pb = max(2, B // 3)
    qb = max(2, B - pb)
    while True:
        p, q = random_prime(pb, rng), random_prime(qb, rng)
        if p != q:
            return p, q


# ---------------------------------------------------------------- log parsing
TS = r"^\d+/\d+/\d+ \d+:\d+:\d+,\s*"
RE_FINISHED = re.compile(TS + r"(?:ecm: )?[Ff]inished (\d+) curves")
RE_COMMENCE = re.compile(TS + r"(?:ecm: )?[Cc]ommencing (\d+) curves")
# also GMP-ECM style: "Finished N curves using GMP-ECM ...", "prp30 = ... (curve 103 stg2 ...)"
# e.g. "prp20 = 7724... (curve=68 stg=2 B1=11000 ...)"  or  "p20 = ..."  or "c30 = ..."
RE_FACTOR = re.compile(TS + r"(?:prp|p|c)(\d+) = (\d+)(?:\s*\(curve[= ](\d+))?")


def parse_segment(text, n):
    """Parse the log text produced by one yafu run on input n."""
    r = dict(found=False, factors=[], curves_run=0, curves_started=0, find_curve=None)
    for line in text.splitlines():
        m = RE_FINISHED.match(line)
        if m:
            r["curves_run"] += int(m.group(1))
            continue
        m = RE_COMMENCE.match(line)
        if m:
            r["curves_started"] += int(m.group(1))
            continue
        m = RE_FACTOR.match(line)
        if m:
            f = int(m.group(2))
            if 1 < f < n and n % f == 0:        # nontrivial divisor of the input
                r["found"] = True
                r["factors"].append(f)
                if m.group(3) and r["find_curve"] is None:
                    r["find_curve"] = int(m.group(3))
    return r


# ------------------------------------------------------------------- driver
def run_one(args, n, logpath, cwd):
    size0 = os.path.getsize(logpath) if os.path.exists(logpath) else 0
    cmd = [args.yafu, f"ecm({n}, {args.C})", "-B1ecm", str(args.B1)] + args.yafu_args
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                              timeout=args.timeout)
        out, rc = proc.stdout, proc.returncode
    except subprocess.TimeoutExpired:
        out, rc = "", -1
        print("  timeout!", file=sys.stderr)
    dt = time.time() - t0
    seg = ""
    if os.path.exists(logpath):
        with open(logpath, "rb") as f:
            f.seek(size0)
            seg = f.read().decode("utf-8", errors="replace")
    res = parse_segment(seg, n)
    res.update(time=dt, rc=rc, logtext=seg, stdout=out)
    return res


def stat_line(name, vals):
    if not vals:
        return f"  {name:<28} n/a"
    return (f"  {name:<28} min={min(vals)}  max={max(vals)}  "
            f"mean={statistics.mean(vals):.1f}  median={statistics.median(vals):.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-B", type=int, required=True, help="bit size of each input")
    ap.add_argument("-N", type=int, required=True, help="number of inputs to test")
    ap.add_argument("-C", type=int, required=True, help="curves per input")
    ap.add_argument("--B1", type=int, required=True, help="ECM stage 1 bound")
    ap.add_argument("--yafu", default="./yafu", help="path to yafu executable")
    ap.add_argument("--log", default=None,
                    help="path to factor.log (default: factor.log in yafu's directory)")
    ap.add_argument("--cwd", default=None,
                    help="working directory to run yafu in (default: yafu's directory)")
    ap.add_argument("--seed", type=int, default=None, help="RNG seed (reproducible inputs)")
    ap.add_argument("--timeout", type=float, default=None, help="per-run timeout, seconds")
    ap.add_argument("--csv", default=None, help="write per-input results to this CSV")
    ap.add_argument("--extra-args", default="",
                    help='string of extra args appended to the yafu call, '
                         'e.g. --extra-args "-threads 4 -v"')
    args = ap.parse_args()
    args.yafu_args = shlex.split(args.extra_args)

    if args.B < 6:
        ap.error("B must be at least 6")
    args.yafu = os.path.abspath(args.yafu)
    cwd = args.cwd or os.path.dirname(args.yafu)
    logpath = args.log or os.path.join(cwd, "factor.log")
    rng = random.Random(args.seed)

    print(f"B={args.B} bits (p~{max(2, args.B // 3)}b, q~{args.B - max(2, args.B // 3)}b)  "
          f"N={args.N}  C={args.C}  B1={args.B1}")
    print(f"yafu={args.yafu}  log={logpath}\n")

    results = []
    for i in range(1, args.N + 1):
        p, q = make_input(args.B, rng)
        n = p * q
        r = run_one(args, n, logpath, cwd)
        r.update(idx=i, n=n, p=p, q=q)
        # which of the known primes did we find?
        r["got_p"] = p in r["factors"]
        r["got_q"] = q in r["factors"]
        results.append(r)
        status = "FOUND" if r["found"] else "no factor"
        print(f"[{i}/{args.N}] {n.bit_length()}b  {status:<9} "
              f"curves={r['curves_run']}  find_curve={r['find_curve']}  "
              f"{r['time']:.2f}s")
        if r["curves_run"] == 0 and not r["found"]:
            print("  warning: no ecm curve lines parsed from log for this run "
                  f"(rc={r['rc']}); check yafu path/log path/log format", file=sys.stderr)

    # ------------------------------------------------------------- summary
    found = [r for r in results if r["found"]]
    nf = len(found)
    total_curves = sum(r["curves_run"] for r in results)
    print("\n================ SUMMARY ================")
    print(f"inputs tested:            {len(results)}  ({args.B}-bit, p={max(2, args.B // 3)}b x q={args.B - max(2, args.B // 3)}b)")
    print(f"factors found:            {nf}  ({100.0 * nf / len(results):.1f}%)")
    print(f"  small prime (B/3) found:  {sum(r['got_p'] for r in results)}")
    print(f"  large prime (2B/3) found: {sum(r['got_q'] for r in results)}")
    print(f"not found:                {len(results) - nf}")
    print(f"total curves run:         {total_curves}  (as reported by yafu; AVX-ECM rounds up to vector size)")
    print(f"avg curves per input:     {total_curves / len(results):.1f}")
    print(stat_line("curves run, all inputs", [r["curves_run"] for r in results]))
    print(stat_line("curves run, found only", [r["curves_run"] for r in found]))
    print(stat_line("curve # of find (curve=)", [r["find_curve"] for r in found
                                                  if r["find_curve"] is not None]))
    if nf:
        print(f"curves per factor found:  {total_curves / nf:.1f}  (total curves / factors found)")
    tt = sum(r["time"] for r in results)
    print(f"total wall time:          {tt:.1f}s  ({tt / len(results):.2f}s per input)")

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("idx,bits,found,curves_run,curves_started,find_curve,got_small,got_large,time,p,q\n")
            for r in results:
                f.write(f"{r['idx']},{r['n'].bit_length()},{int(r['found'])},{r['curves_run']},"
                        f"{r['curves_started']},{'' if r['find_curve'] is None else r['find_curve']},"
                        f"{int(r['got_p'])},{int(r['got_q'])},{r['time']:.3f},{r['p']},{r['q']}\n")
        print(f"per-input results written to {args.csv}")


if __name__ == "__main__":
    main()
