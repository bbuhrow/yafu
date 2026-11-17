#!/usr/bin/env python3
"""Fit whole-job GNFS times T(D) for -plan bayesian and emit "nfs <digits> <s>"
lines for bayes_tune.txt (docs/ecm_bayes.md).

Usage:
  python fit_sieve_times.py RUNS.csv [--testsieve TS.csv] [--threads N]
                            [--from 110] [--to 200] [--step 5] [--model L|linear]
                            [--append FILE]

RUNS.csv   output of parse_yafu_nfs.py: complete local yafu NFS runs (digits,
           total_s, sieve_threads, poly_s, sieve_filter_s, la_s, sqrt_s,
           matrix_cols ...).  Two or more sizes are needed for a fit.
TS.csv     optional test-sieve results for larger sizes, columns
           digits,rels_needed,rels_per_sec   (rels_per_sec at the same thread
           count as the runs, averaged over the q range; rels_needed = the raw
           relations yafu would gather, e.g. its min_rels).  For these sizes
           T = poly + rels_needed / rels_per_sec + post, with poly and post
           (filtering + LA + sqrt) extrapolated from the runs.

Model (fitted in log space, least squares):
  complete runs:     ln T = a + b D
  post-processing:   ln post = c + d D   (LA ~ matrix^1.75, matrix ~1.07^D
                     from the NFS@Home medians C150-C175; with one run only,
                     d = 1.75 ln 1.07 = 0.118 is assumed)
  poly selection:    fraction of T from the runs (yafu's poly budget scales
                     with the expected sieve time).  It is used as logged:
                     yafu logs "poly selection with 1 threads" (msieve threads
                     internally) and stops at a wall-clock deadline.
When test-sieve rows exist they replace the run fit at their sizes, and the
emitted curve is fitted over runs + test-sieve points together:
  --model L (default):  ln T = a + k * L(D), L = (64/9)^(1/3) x^(1/3) (ln x)^(2/3),
                        x = D ln 10: the GNFS L[1/3] shape with a fitted scale k,
                        so the doubling distance grows with size (C120-C140 runs
                        on an i5-11600K: x2 per 6.3 digits at C140, 7.3 at C180)
  --model linear:       ln T = a + b D (x2 per 6.0 digits everywhere; this
                        overshoots the L shape by ~1.7x at C180)
Both fit C120-C140 equally well; they differ only in extrapolation, which
test-sieve points at C150+ pin down.

Each emitted line is "nfs <digits> <seconds> <threads>"; yafu scales the time
linearly when run with a different thread count.
"""
import argparse
import csv
import math
import sys


def linfit(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        raise ValueError('need at least two different sizes')
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return my - b * mx, b


def lshape(D):
    x = D * math.log(10.0)
    return (64.0 / 9.0) ** (1.0 / 3.0) * x ** (1.0 / 3.0) * math.log(x) ** (2.0 / 3.0)


def num(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs')
    ap.add_argument('--testsieve')
    ap.add_argument('--threads', type=int, default=0,
                    help='thread count the times are for (default: the runs\' sieve threads)')
    ap.add_argument('--from', dest='d0', type=int, default=110,
                    help='first size written (yafu does not extrapolate more than 5 digits '
                    'below it, and uses the tune() fit there instead)')
    ap.add_argument('--to', dest='d1', type=int, default=200)
    ap.add_argument('--step', type=int, default=5)
    ap.add_argument('--model', choices=('L', 'linear'), default='L')
    ap.add_argument('--append', help='append the lines to this bayes_tune.txt '
                    '(replacing its existing nfs lines)')
    a = ap.parse_args()

    runs = [r for r in csv.DictReader(open(a.runs)) if num(r, 'total_s') and num(r, 'digits')]
    if not runs:
        sys.exit('no complete runs in ' + a.runs)
    thr = {int(num(r, 'sieve_threads') or 0) for r in runs}
    if len(thr) > 1:
        print('warning: runs use different thread counts %s' % sorted(thr), file=sys.stderr)
    threads = a.threads or max(thr)

    pts = [(num(r, 'digits'), num(r, 'total_s')) for r in runs]
    post = []
    polyfrac = []
    for r in runs:
        p = sum(v for v in (num(r, 'la_s'), num(r, 'sqrt_s'), num(r, 'relproc_s')) if v)
        if p:
            post.append((num(r, 'digits'), p))
        if num(r, 'poly_s'):
            polyfrac.append(num(r, 'poly_s') / num(r, 'total_s'))
    pf = sum(polyfrac) / len(polyfrac) if polyfrac else 0.08

    if len(post) >= 2:
        c, d = linfit([x for x, _ in post], [math.log(y) for _, y in post])
    elif post:
        d = 1.75 * math.log(1.07)
        c = math.log(post[0][1]) - d * post[0][0]
    else:
        c = d = None

    if a.testsieve:
        for r in csv.DictReader(open(a.testsieve)):
            D, need, rate = num(r, 'digits'), num(r, 'rels_needed'), num(r, 'rels_per_sec')
            if not (D and need and rate):
                continue
            sieve = need / rate
            pp = math.exp(c + d * D) if c is not None else 0.25 * sieve
            T = (sieve + pp) / (1.0 - pf)
            pts.append((D, T))
            print('# test-sieve C%.0f: sieve %.4g s + post %.4g s, poly %.0f%% -> %.4g s'
                  % (D, sieve, pp, 100 * pf, T), file=sys.stderr)

    if len({x for x, _ in pts}) < 2:
        sys.exit('need runs/test-sieve points at two or more sizes')
    f = lshape if a.model == 'L' else (lambda D: D)
    A, B = linfit([f(x) for x, _ in pts], [math.log(y) for _, y in pts])
    model = lambda D: A + B * f(D)
    resid = [math.log(y) - model(x) for x, y in pts]
    sd = math.sqrt(sum(e * e for e in resid) / max(1, len(resid) - 2)) if len(resid) > 2 else 0.0
    dbl = lambda D: math.log(2) / (model(D + 0.5) - model(D - 0.5))

    lo, hi = min(x for x, _ in pts), max(x for x, _ in pts)
    lines = ['# sieve times: fit_sieve_times.py, %s model, %d points C%.0f-C%.0f, '
             'x2 per %.1f digits at C%.0f / %.1f at C%d, resid sd %.2f (ln), %d threads'
             % (a.model, len(pts), lo, hi, dbl(hi), hi, dbl(a.d1), a.d1, sd, threads)]
    for D in range(a.d0, a.d1 + 1, a.step):
        lines.append('nfs %d %.6g %d' % (D, math.exp(model(D)), threads))
    out = '\n'.join(lines) + '\n'

    if a.append:
        try:
            old = open(a.append).read().splitlines(True)
        except FileNotFoundError:
            old = []
        keep = [l for l in old if not (l.startswith('nfs ') or l.startswith('# sieve times'))]
        hdr_thr = [l for l in keep if l.startswith('threads:')]
        if hdr_thr and int(hdr_thr[0].split(':')[1]) != threads:
            print('warning: %s is for %s threads, these times for %d'
                  % (a.append, hdr_thr[0].split(':')[1].strip(), threads), file=sys.stderr)
        open(a.append, 'w').write(''.join(keep) + out)
        print('wrote %d nfs lines to %s' % (len(lines) - 1, a.append), file=sys.stderr)
    else:
        sys.stdout.write(out)


if __name__ == '__main__':
    main()
