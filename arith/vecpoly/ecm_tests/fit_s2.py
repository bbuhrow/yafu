"""
Fit the FFT stage-2 cost model (avxecm_stage2fft.c, s2_cost) to
vp_test -b output, per ISA.

  python3 fit_s2.py bench_ifma.txt bench_f.txt bench_avx2.txt

Model (units: template mulmods, mu = one 8-lane template mulmod):
  T  = tc(t, n) dF,  tc = a t + b + c t^2 / n     G's product tree
  per block  : 10 dF (giant roots) + T (1 + mm)   mm = mulmod / tree
  per batch  : 17 dF + T (1 + mi + ev + rebuilt)  mi = 1/F + transforms,
                                                  ev = evaluation
Lean mode (no stored transforms of 1/F, F; the benchmark's lean also keeps
only every 3rd tree level) gives mm_lean, mi_lean and ev_lean (= lean
evaluation minus the stride-3 rebuild).
"""
import re
import sys
import numpy as np


def parse(fn):
    rows = []
    n = None
    for line in open(fn):
        m = re.match(r'(\d+) bits: n = (\d+)', line)
        if m:
            n = int(m.group(2))
            continue
        m = re.search(r'(full|lean) m=\s*(\d+): tree\s+([\d.]+) ms, roots->poly\s+([\d.]+) ms, '
                      r'modf_init\s+([\d.]+) ms, mulmod\s+([\d.]+) ms, eval\s+([\d.]+) ms.*'
                      r'mulmod = (\d+) ns; tree = ([\d.]+)', line)
        if m:
            M = int(m.group(2))
            t = M.bit_length() - 1
            tree, rp, mi, mm, ev, mu, tc = [float(x) for x in m.groups()[2:]]
            rows.append(dict(kind=m.group(1), n=n, t=t, tree=tree, rp=rp, mi=mi, mm=mm,
                             ev=ev, mu=mu, tc=tc))
    return rows


def rebuilt(t, k):
    return sum(s % k for s in range(1, t)) / t


def fit(fn):
    rows = parse(fn)
    full = [r for r in rows if r['kind'] == 'full' and r['t'] >= 8]
    A = np.array([[r['t'], 1.0, r['t'] ** 2 / r['n']] for r in full])
    y = np.array([r['tc'] for r in full])
    (a, b, c), _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    err = max(abs(A[i] @ [a, b, c] / y[i] - 1) for i in range(len(y)))
    big = [r for r in rows if r['t'] >= 9]
    med = lambda kind, key: float(np.median([r[key] / r['rp'] for r in big if r['kind'] == kind]))
    ev_lean = float(np.median([r['ev'] / r['rp'] - rebuilt(r['t'], 3)
                               for r in big if r['kind'] == 'lean']))
    mus = {}
    for r in rows:
        mus.setdefault(r['n'], []).append(r['mu'])
    mus = {k: float(np.median(v)) for k, v in mus.items()}
    return dict(a=a, b=b, c=c, maxerr=err, mm=med('full', 'mm'), mm_lean=med('lean', 'mm'),
                mi=med('full', 'mi'), mi_lean=med('lean', 'mi'), ev=med('full', 'ev'),
                ev_lean=ev_lean, mu=mus)


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        f = fit(fn)
        print('%s: tc = %.2f t %+.1f %+.2f t^2/n (max err %.0f%%)' % (fn, f['a'], f['b'], f['c'], 100 * f['maxerr']))
        print('   mm %.2f / lean %.2f, mi %.2f / lean %.2f, ev %.2f / lean %.2f'
              % (f['mm'], f['mm_lean'], f['mi'], f['mi_lean'], f['ev'], f['ev_lean']))
        print('   mu(n) ns: ' + ', '.join('%d: %.0f' % kv for kv in sorted(f['mu'].items())))
        print('   C: { %.3f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f },'
              % (f['a'], f['b'], f['c'], f['mm'], f['mm_lean'], f['mi'], f['mi_lean'], f['ev'], f['ev_lean']))
