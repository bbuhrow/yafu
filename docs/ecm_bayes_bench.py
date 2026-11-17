#!/usr/bin/env python3
"""Benchmark factor() pretest plans on sets of random composites
(docs/ecm_bayes.md: normal vs bayesian).

Usage:
  python ecm_bayes_bench.py --yafu ./yafu [--sizes 70,80,90] [--count 20]
      [--plans normal,bayesian] [--threads 6] [--kind random|semiprime]
      [--seed 1] [--workdir DIR] [--out results.csv] [-- extra yafu args]

For each size it draws --count composites (random D-digit composites, or
balanced semiprimes with --kind semiprime), factors the same numbers with each
plan (one yafu process per plan and size, numbers fed as a batch on stdin), and
reports per-plan total and median time plus the paired difference with its
standard error.  Run it on an otherwise idle machine: the comparison is only
as good as the timing noise.  The work directory needs yafu.ini (with tune_info)
and bayes_tune.txt for the bayesian plan; it defaults to yafu's own directory.
"""
import argparse
import csv
import math
import os
import random
import re
import statistics
import subprocess
import sys
import time


def is_probable_prime(n, rounds=20):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for _ in range(rounds):
        a = random.randrange(2, n - 1)
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


def rand_prime(digits):
    while True:
        p = random.randrange(10 ** (digits - 1), 10 ** digits) | 1
        if is_probable_prime(p):
            return p


def make_inputs(kind, digits, count):
    out = []
    while len(out) < count:
        if kind == 'semiprime':
            a = digits // 2
            n = rand_prime(a) * rand_prime(digits - a)
            if len(str(n)) != digits:
                continue
        else:
            n = random.randrange(10 ** (digits - 1), 10 ** digits)
            if is_probable_prime(n):
                continue
        out.append(n)
    return out


def run_plan(yafu, workdir, plan, threads, nums, extra):
    for f in ('siqs.dat', 'factor.log'):
        try:
            os.remove(os.path.join(workdir, f))
        except FileNotFoundError:
            pass
    inp = ''.join('factor(%d)\n' % n for n in nums)
    t0 = time.time()
    p = subprocess.run([yafu, '-threads', str(threads), '-plan', plan] + extra,
                       input=inp, capture_output=True, text=True, cwd=workdir)
    wall = time.time() - t0
    times = [float(x) for x in re.findall(r'Total factoring time = ([\d.]+) seconds', p.stdout)]
    if len(times) != len(nums):
        sys.exit('plan %s: got %d timings for %d inputs; yafu output tail:\n%s'
                 % (plan, len(times), len(nums), p.stdout[-2000:]))
    return times, wall


def main():
    argv = sys.argv[1:]
    extra = []
    if '--' in argv:
        i = argv.index('--')
        argv, extra = argv[:i], argv[i + 1:]
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--yafu', required=True)
    ap.add_argument('--sizes', default='70,80,90')
    ap.add_argument('--count', type=int, default=20)
    ap.add_argument('--plans', default='normal,bayesian')
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--kind', choices=('random', 'semiprime'), default='random')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--workdir')
    ap.add_argument('--out', default='ecm_bayes_bench.csv')
    a = ap.parse_args(argv)

    yafu = os.path.abspath(a.yafu)
    workdir = a.workdir or os.path.dirname(yafu)
    plans = a.plans.split(',')
    random.seed(a.seed)

    rows = []
    for D in [int(x) for x in a.sizes.split(',')]:
        nums = make_inputs(a.kind, D, a.count)
        res = {}
        for plan in plans:
            res[plan], wall = run_plan(yafu, workdir, plan, a.threads, nums, extra)
            print('C%d %-9s total %8.1f s  median %7.2f s  (wall %.1f s)'
                  % (D, plan, sum(res[plan]), statistics.median(res[plan]), wall), flush=True)
        for i, n in enumerate(nums):
            rows.append(dict(digits=D, n=n, **{p: res[p][i] for p in plans}))
        if len(plans) == 2:
            d = [y - x for x, y in zip(res[plans[0]], res[plans[1]])]
            se = statistics.stdev(d) / math.sqrt(len(d)) if len(d) > 1 else float('nan')
            print('C%d %s - %s: mean %+.2f s per input (se %.2f), total %+.1f%%'
                  % (D, plans[1], plans[0], statistics.mean(d), se,
                     100 * (sum(res[plans[1]]) / sum(res[plans[0]]) - 1)), flush=True)

    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['digits', 'n'] + plans)
        w.writeheader()
        w.writerows(rows)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
