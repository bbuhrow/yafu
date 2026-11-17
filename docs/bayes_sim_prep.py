#!/usr/bin/env python3
# ========== Note: created by Claude with prompting by bbuhrow ==========
"""Inputs and cost model for yafu's bayessim() replay (docs/ecm_bayes.md s. 6).

  prep:     bayes_sim_prep.py prep factor.json [more.json ...] [--out-dir D]
            writes bayes_sim_in.txt (each distinct input, all its prime
            factors) and bayes_sim_cost.txt (the machine's costs, fitted
            from the runs in the json files).
  check:    bayes_sim_prep.py check factor.json --sim-dir D
            compares the replay with the runs it was fitted from: the ECM
            steps of the normal plan, and the mean times per size.

Cost model (all wall clock, as the benchmark ran: same build, threads):
  ECM call     ecm_step + curves * c (B1/1e5)^e (digits/70)^q
               (grid search on e and q, least squares for the rest)
  P-1          a + b digits per level (B1 = 150k, 3.75M, 15M)
  sieve        median SIQS time per composite size, smoothed (log-linear
               local fit), from runs where the SIQS input is known: exactly
               two prime factors above the rho bound and SIQS ran
Runs with a stall (a phase taking > 5 s more than the model) are left out.
"""
import argparse
import collections
import json
import math
import os
import statistics
import sys

import numpy as np

PM1_B1 = (150000, 3750000, 15000000)


def load(paths):
    recs = []
    for p in paths:
        with open(p) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    recs.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    # the json "siqs" runtime is not reset between inputs: SIQS ran on an
    # input when the value changed from the previous record of the same run
    # (same argument string)
    prev = {}
    for r in recs:
        key = r.get('input-argument-string', '')
        s = r.get('runtime', {}).get('siqs')
        r['_siqs'] = s if (s is not None and s != prev.get(key)) else 0.0
        prev[key] = s
    return recs


def plan_of(r):
    a = r.get('input-argument-string', '')
    return 'bayes' if 'bayes' in a else 'normal'


def ecm_digits(r, rho_bound):
    big = [int(p) for p in r['factors-prime'] if int(p) >= rho_bound]
    return len(str(math.prod(big))) if big else 0


def fit_ecm(recs, rho_bound):
    data = []
    for r in recs:
        rt = r['runtime']
        if 'ecm-curves' not in r:
            continue
        data.append(([(int(b), n) for b, n in r['ecm-curves'].items()],
                     ecm_digits(r, rho_bound), rt.get('ecm', 0.0)))

    def solve(rows):
        y = np.array([t for _, _, t in rows])
        best = None
        for e in np.arange(0.4, 1.31, 0.02):
            for q in np.arange(0.0, 3.01, 0.1):
                X = np.array([[len(c), sum(n * (b / 1e5) ** e for b, n in c) * (d / 70) ** q]
                              for c, d, _ in rows])
                co = np.linalg.lstsq(X, y, rcond=None)[0]
                sse = float(np.sum((X @ co - y) ** 2))
                if best is None or sse < best[0]:
                    best = (sse, e, q, co)
        sse, e, q, co = best
        r2 = 1 - sse / float(np.sum((y - y.mean()) ** 2))
        return dict(step=max(co[0], 0.0), c=co[1], e=e, q=q, r2=r2)

    m = solve(data)

    def pred(c, d):
        return m['step'] * len(c) + sum(n * m['c'] * (b / 1e5) ** m['e'] * (d / 70) ** m['q']
                                        for b, n in c)
    clean = [x for x in data if x[2] - pred(x[0], x[1]) < 5.0]
    if len(clean) < len(data):
        m = solve(clean)
    m['n'] = len(clean)
    m['b1max'] = max(b for c, _, _ in clean for b, n in c if n > 0)
    m['dropped'] = len(data) - len(clean)
    return m


def fit_pm1(recs, rho_bound):
    rows, y = [], []
    for r in recs:
        if 'pm1-curves' not in r:
            continue
        d = ecm_digits(r, rho_bound)
        x = []
        for b in PM1_B1:
            n = 1 if str(b) in r['pm1-curves'] else 0
            x += [n, n * d]
        if r['runtime'].get('pm1', 0) > 5 + 2 * sum(x[2:]):     # stall
            continue
        rows.append(x)
        y.append(r['runtime'].get('pm1', 0.0))
    A = np.array(rows, float)
    y = np.array(y)
    out = {}
    used = [j for j in range(A.shape[1]) if A[:, j].any()]
    co = np.linalg.lstsq(A[:, used], y, rcond=None)[0]
    full = np.zeros(A.shape[1])
    full[used] = co
    for i, b in enumerate(PM1_B1):
        if A[:, 2 * i].any():
            out[b] = (full[2 * i], full[2 * i + 1])
    # levels never reached: scale the highest measured level by B1
    for i, b in enumerate(PM1_B1):
        if b not in out:
            lower = [bb for bb in out if bb < b]
            if lower:
                bb = max(lower)
                s = b / bb
                out[b] = (out[bb][0] * s, out[bb][1] * s)
    return out


def fit_siqs(recs, rho_bound):
    by = collections.defaultdict(list)
    for r in recs:
        if not r['_siqs']:
            continue
        big = [int(p) for p in r['factors-prime'] if int(p) >= rho_bound]
        if len(big) != 2:
            continue
        by[len(str(big[0] * big[1]))].append(r['_siqs'])
    pts = sorted((d, statistics.median(v), len(v)) for d, v in by.items() if len(v) >= 3)
    if len(pts) < 3:
        return pts, []
    ds = np.array([p[0] for p in pts], float)
    lt = np.log([p[1] for p in pts])
    w = np.sqrt([p[2] for p in pts])
    sm = []
    for d in ds:
        k = np.abs(ds - d) <= 3
        co = np.polyfit(ds[k] - d, lt[k], 1, w=w[k]) if k.sum() >= 3 else (0, lt[ds == d][0])
        sm.append((int(d), float(math.exp(co[-1]))))
    # extend 10 digits each way with the slope of the outer 8 digits
    ext = []
    for side in (0, 1):
        seg = sm[:8] if side == 0 else sm[-8:]
        sl = np.polyfit([a for a, _ in seg], [math.log(b) for _, b in seg], 1)[0]
        d0, t0 = sm[0] if side == 0 else sm[-1]
        for k in (5, 10):
            dd = d0 - k if side == 0 else d0 + k
            ext.append((dd, t0 * math.exp(sl * (dd - d0))))
    sm = sorted(sm + ext)
    return pts, sm


def cmd_prep(a):
    recs = load(a.json)
    rho = a.rho_bound
    inputs = {}
    for r in recs:
        n = r['input-decimal']
        if n not in inputs:
            inputs[n] = r['factors-prime']
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, 'bayes_sim_in.txt'), 'w') as f:
        f.write('# <class> <prime factors> ... from %s\n' % ' '.join(a.json))
        for n, fs in inputs.items():
            prod = math.prod(int(p) for p in fs)
            if prod != int(n):
                print('skipping input whose factors do not multiply out:', n[:30])
                continue
            f.write('%d %s\n' % (len(n), ' '.join(fs)))

    ecm = fit_ecm(recs, rho)
    pm1 = fit_pm1(recs, rho)
    pts, sq = fit_siqs(recs, rho)
    print('ecm: step %.4f s + curves * %.3e (B1/1e5)^%.2f (d/70)^%.1f  R2 %.3f (%d runs, %d stalls dropped)'
          % (ecm['step'], ecm['c'], ecm['e'], ecm['q'], ecm['r2'], ecm['n'], ecm['dropped']))
    for b, (x, y) in sorted(pm1.items()):
        print('pm1 B1=%d: %.4f + %.5f d s' % (b, x, y))
    print('siqs: %d sizes with data, C%d-C%d' % (len(pts), pts[0][0], pts[-1][0]) if pts else 'siqs: no data')

    with open(os.path.join(a.out_dir, 'bayes_sim_cost.txt'), 'w') as f:
        f.write('# bayessim cost model fitted by bayes_sim_prep.py from %s\n' % ' '.join(a.json))
        f.write('reps %d\nseed %d\nlanes %d\nrho_bound %g\n' % (a.reps, a.seed, a.lanes, rho))
        f.write('# ECM: per call + per curve c (B1/1e5)^e (digits/70)^q; R2 %.3f over %d runs\n'
                % (ecm['r2'], ecm['n']))
        f.write('ecm_step %.6g\necm_curve %.6g %.4g %.4g\n' % (ecm['step'], ecm['c'], ecm['e'], ecm['q']))
        f.write('# above the largest fitted B1 the per-curve time grows as B1^ecm_e_hi (stage 1 is linear in B1)\n')
        f.write('ecm_b1max %d\necm_e_hi %.3g\n' % (ecm['b1max'], a.ecm_e_hi))
        f.write('# P-1: B1 B2 a b (seconds a + b digits); B2 approximates GMP-ECM P-1 default\n')
        for b, (x, y) in sorted(pm1.items()):
            b2 = a.pm1_b2_scale * (11.0 / 6.0 * b) ** 1.43
            f.write('pm1 %d %.6g %.6g %.6g\n' % (b, b2, max(x, 0.0), y))
        f.write('# sieve: digits seconds (smoothed medians; first/last two are extrapolated)\n')
        for d, t in sq:
            f.write('siqs %d %.6g\n' % (d, t))
        f.write('# plans: the first is the baseline of the comparison\n')
        for p in a.policy:
            f.write('policy %s\n' % p)
    print('wrote %s/bayes_sim_in.txt (%d inputs) and bayes_sim_cost.txt' % (a.out_dir, len(inputs)))


def parse_path(s):
    out = collections.Counter()
    for tok in s.split(';'):
        if ':' in tok:
            b, n = tok.split(':')
            out[int(b)] += int(n)
    return out


def cmd_check(a):
    recs = load(a.json)
    sim_idx = {}
    with open(os.path.join(a.sim_dir, 'bayes_sim_in.txt')) as f:
        i = 0
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            sim_idx[str(math.prod(int(p) for p in line.split()[1:]))] = i
            i += 1
    paths = {}
    with open(os.path.join(a.sim_dir, 'bayes_sim_paths.csv')) as f:
        next(f)
        for line in f:
            plan, cls, idx, path = line.rstrip('\n').split(',', 3)
            paths[(plan, int(idx))] = parse_path(path)
    sims = collections.defaultdict(dict)
    with open(os.path.join(a.sim_dir, 'bayes_sim_out.csv')) as f:
        next(f)
        for line in f:
            plan, cls, idx, ecm, pm1, sieve, total, tl = line.strip().split(',')
            sims[plan][int(idx)] = (int(cls), float(ecm), float(pm1), float(sieve), float(total))

    # 1. ECM steps of the normal plan: the observed curves are the all-fail
    # path, or a prefix of it when a step found a factor
    stats = collections.defaultdict(lambda: [0, 0, 0])
    for r in recs:
        if plan_of(r) != 'normal' or 'ecm-curves' not in r:
            continue
        i = sim_idx.get(r['input-decimal'])
        if i is None or ('normal', i) not in paths:
            continue
        obs = collections.Counter({int(b): n for b, n in r['ecm-curves'].items()})
        sim = paths[('normal', i)]
        cls = len(r['input-decimal'])
        st = stats[cls]
        if obs == sim:
            st[0] += 1
        elif all(obs[b] <= sim[b] for b in obs):
            st[1] += 1
        else:
            st[2] += 1
            if a.verbose:
                print('  C%d mismatch: observed %s simulated %s' % (cls, dict(obs), dict(sim)))
    print('normal plan ECM steps, observed vs replay (all-fail path):')
    for cls in sorted(stats):
        e, p, m = stats[cls]
        print('  C%d: identical %d, observed shorter (a find) %d, different %d' % (cls, e, p, m))

    # 2. mean times per class: observed (stall runs dropped) vs replay
    print('mean seconds per input, observed (ecm + pm1 + siqs) vs replay:')
    for plan_sim, plan_obs in (('normal', 'normal'), (a.bayes_plan, 'bayes')):
        if plan_sim not in sims:
            continue
        obs = collections.defaultdict(list)
        for r in recs:
            if plan_of(r) != plan_obs:
                continue
            rt = r['runtime']
            if rt['total'] - rt.get('ecm', 0) - rt.get('pm1', 0) - r['_siqs'] > 5 or rt.get('ecm', 0) > 5:
                continue
            obs[len(r['input-decimal'])].append((rt.get('ecm', 0), rt.get('pm1', 0), r['_siqs']))
        simc = collections.defaultdict(list)
        for i, v in sims[plan_sim].items():
            simc[v[0]].append(v[1:4])
        for cls in sorted(obs):
            o = np.mean(np.array(obs[cls]), axis=0)
            s = np.mean(np.array(simc[cls]), axis=0) if simc[cls] else np.zeros(3)
            print('  %-8s C%d  observed ecm %.4f pm1 %.4f siqs %.4f = %.4f | replay ecm %.4f pm1 %.4f '
                  'sieve %.4f = %.4f' % (plan_sim, cls, o[0], o[1], o[2], o.sum(), s[0], s[1], s[2], s.sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prep')
    p.add_argument('json', nargs='+')
    p.add_argument('--out-dir', default='.')
    p.add_argument('--rho-bound', type=float, default=1e6)
    p.add_argument('--reps', type=int, default=16)
    p.add_argument('--seed', type=int, default=12345)
    p.add_argument('--lanes', type=int, default=8, help='AVX-ECM curves per vector (8 for AVX-512)')
    p.add_argument('--pm1-b2-scale', type=float, default=3.0)
    p.add_argument('--ecm-e-hi', type=float, default=1.0,
                   help='B1 exponent of the ECM curve time above the largest fitted B1')
    p.add_argument('--policy', action='append',
                   help='plan line: "<name> normal|light|deep|noecm|custom <r>|bayes <siqs_scale> <ecm_scale>|bayes_tune"')
    c = sub.add_parser('check')
    c.add_argument('json', nargs='+')
    c.add_argument('--sim-dir', default='.')
    c.add_argument('--bayes-plan', default='bayes', help='replay plan to set against the observed bayesian runs')
    c.add_argument('-v', '--verbose', action='store_true')
    a = ap.parse_args()
    if a.cmd == 'prep':
        if not a.policy:
            a.policy = ['normal normal', 'bayes bayes 1 1', 'bayes_sq0.5 bayes 0.5 1',
                        'bayes_sq2 bayes 2 1', 'light light', 'deep deep']
        a.rho_bound = int(a.rho_bound)
        cmd_prep(a)
    else:
        cmd_check(a)


if __name__ == '__main__':
    sys.exit(main())
