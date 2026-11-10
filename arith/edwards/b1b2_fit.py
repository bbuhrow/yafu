#!/usr/bin/env python3
"""Refit the per-curve cost model of b1b2_study.py to uedw_lvlbench
measurements: per input, time = INPUT + sum over curves run of
(a + b M1 + [stage 2 run] (c + d X)), X = crosses + 3.5 giants, with the
curve counts (and which curves ran stage 2) from the exact model.  Least
squares on (a, b, c, d) with c fixed in proportion (setup ~ 100 crosses).
usage: b1b2_fit.py b1b2_lvlbench.txt"""
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b1b2_study as st         # noqa: E402
import b1b2_emit as em          # noqa: E402

V2 = {B1: (pol, m) for B1, pol, m in em.LEVELS2}


def counts(rows, s, st2, K=48):
    """mean per input: curves, stage-2 runs, GKL builds (+ starts)"""
    nc = n2 = ng = 0
    for _, ords in rows:
        for k in range(K):
            if k >= st.NTAB:
                ng += 1
            out = st.curve_outcome(ords[k], s, st2)
            nc += 1
            if out == 1:
                break
            n2 += 1
            if out == 2:
                break
    n = len(rows)
    return nc / n, n2 / n, ng / n


def main():
    data = []
    rowsN = {}
    for line in open(sys.argv[1]):
        m = re.match(r'N (\d+) set (\d) B1\s+(\d+):\s+\d+ ns per find\s+(\d+) ns per input', line)
        if not m:
            continue
        N, sset, B1, t = int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4))
        if N not in rowsN:
            rowsN[N] = st.load(N)
        pol, mult = ('lvl', st.legacy_mult(B1)) if sset == 1 else V2[B1]
        s = st.get_scalar(B1, pol)
        s2o = st.stage2(B1, mult, pol)
        M1 = st.chain_m(s)
        X = (s2o.cost - st.S2_FIX) / st.S2_X           # crosses + 3.5 giants
        nc, n2, ng = counts(rowsN[N], s, s2o)
        data.append((t, nc, nc * M1, n2, n2 * X, ng, N, sset, B1))
    # unknowns: a (per curve), b (per M1), c (per stage 2), d (per X); fixed:
    # input 80 ns, GKL build 78 ns
    import itertools
    best = None
    for a in range(0, 401, 10):
        for c in range(0, 601, 10):
            # b, d by least squares given a, c
            sxx = sxy = syy = sx1 = sy1 = 0.0
            A11 = A12 = A22 = r1 = r2 = 0.0
            for t, nc, ncM, n2, n2X, ng, *_ in data:
                y = t - 80 - 78 * ng - a * nc - c * n2
                A11 += ncM * ncM
                A12 += ncM * n2X
                A22 += n2X * n2X
                r1 += ncM * y
                r2 += n2X * y
            det = A11 * A22 - A12 * A12
            b = (r1 * A22 - r2 * A12) / det
            d = (A11 * r2 - A12 * r1) / det
            err = 0.0
            for t, nc, ncM, n2, n2X, ng, *_ in data:
                pred = 80 + 78 * ng + a * nc + b * ncM + c * n2 + d * n2X
                err += (math.log(pred / t)) ** 2
            if best is None or err < best[0]:
                best = (err, a, b, c, d)
    err, a, b, c, d = best
    print('fit: per curve %.0f + %.3f ns/M1; per stage 2 %.0f + %.3f ns/X; rms log err %.3f'
          % (a, b, c, d, math.sqrt(err / len(data))))
    for t, nc, ncM, n2, n2X, ng, N, sset, B1 in data:
        pred = 80 + 78 * ng + a * nc + b * ncM + c * n2 + d * n2X
        print('  N %d set %d B1 %3d: real %7.0f fit %7.0f  %.3f' % (N, sset, B1, t, pred, t / pred))


if __name__ == '__main__':
    main()
