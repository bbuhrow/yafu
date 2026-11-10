#!/usr/bin/env python3
"""
Exact oracle for microecm_edw's re-tuned parameter set 2 (edw_tiny_v2.h,
plan doc Part I): one curve per (n = p q, curve id) at every v2 level, the
stage of the find vs the prediction from exact ord P mod p and mod q, the
level's s and an exact replay of the generic stage-2 walk (giant windows
A = 30 + 120 j from [30]Q, cross terms per table entry, stop before a
product 0 mod n).

  uedw_oracle2.py curves BIN [nper] [seed]
"""
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uedw_oracle as uo        # noqa: E402
import b1b2_study as st         # noqa: E402
import b1b2_emit as em          # noqa: E402
import bs_oracle as bo          # noqa: E402
import edw_ref as er            # noqa: E402

MAP = uo.MAP


def gen_zero_flags(P0, A, p, s2):
    """per table cross term, whether it is 0 mod p (the C walk exactly)"""
    a24 = (A + 2) * pow(4, -1, p) % p

    def dup(P):
        X, Z = P
        U = (X - Z) ** 2 % p
        V = (X + Z) ** 2 % p
        w = (V - U) % p
        return (U * V % p, w * (U + a24 * w) % p)

    def uadd(P1, P2, Pin):
        u = (P1[0] - P1[1]) * (P2[0] + P2[1]) % p
        v = (P1[0] + P1[1]) * (P2[0] - P2[1]) % p
        return (Pin[1] * (u + v) ** 2 % p, Pin[0] * (u - v) ** 2 % p)

    Pb = [None] * 18
    Pbprod = [0] * 18
    Pb[1] = P0
    Pbprod[1] = P0[0] * P0[1] % p
    Pb[2] = dup(Pb[1])
    Pb[3] = uadd(Pb[1], Pb[2], Pb[1])
    pt6 = dup(Pb[3])
    pt5 = uadd(Pb[3], Pb[2], Pb[1])
    Pb[3] = pt5
    Pb[4] = uadd(pt6, pt5, Pb[1])
    for (i, k, j) in ((3, 4, 5), (1, 3, 1)):
        if i == 1:
            Pb[3] = uadd(pt6, Pb[1], pt5)
        while j + 12 < 60:
            Pb[MAP[j + 12]] = uadd(pt6, Pb[k], Pb[i])
            i, k, j = k, MAP[j + 12], j + 6
    Pd = uadd(Pb[9], Pb[10], Pb[2])
    pt4 = dup(Pb[2])
    for t in range(3, 18):
        Pbprod[t] = Pb[t][0] * Pb[t][1] % p
    Pad = uadd(Pb[MAP[17]], Pb[MAP[13]], pt4)
    Pa = uadd(Pd, Pad, Pad)
    Pd = uadd(Pa, Pad, Pd)
    cur, prev = Pad, Pa
    flags = []
    for b in s2:
        if b == 0:
            cur, prev = uadd(cur, Pd, prev), cur
            continue
        Rprod = cur[0] * cur[1] % p
        tt1 = (cur[0] - Pb[MAP[b]][0]) * (cur[1] + Pb[MAP[b]][1]) % p
        flags.append((tt1 + Pbprod[MAP[b]] - Rprod) % p == 0)
    return flags


def predict(p, cid, s, s2):
    try:
        d, x, y = uo.curve_ref(cid, p)
        if (1 + d) % p == 0:
            raise er.Find()
        A, g, u, v = er.edw_to_mont(d, x, y, p)
    except (er.Find, ValueError, ZeroDivisionError):
        return None
    R0, _ = bo.x_ladder(s, u, A, p)
    return (R0[1] % p == 0, ([], gen_zero_flags(R0, A, p, s2)))


# factor size per level (the v2 dispatch: balanced n)
LEVEL_BITS = {20: 17, 27: 19, 47: 21, 70: 24, 85: 26, 125: 28, 165: 30, 205: 32}


def check(binp, nper=200, seed=1):
    d = uo.Driver(binp)
    d.ask('v 2')
    rng = random.Random(seed)
    tot = bad = degen = exc = 0
    stats = {}
    for B1, pol, mult in em.LEVELS2:
        s = st.scalar(B1, pol)
        s2 = em.gen_table(B1, mult * B1)
        pb = LEVEL_BITS[B1]
        for _ in range(nper):
            p = uo.rand_prime(rng, pb)
            q = uo.rand_prime(rng, pb + rng.randrange(0, 4))
            if p == q or p * q >= 1 << 64:
                continue
            n = p * q
            m = rng.randrange(1, 400)
            if rng.random() < 0.3:
                m = uo.TABLE | rng.randrange(len(uo.table_curves()))
            sp, sq = predict(p, m, s, s2), predict(q, m, s, s2)
            stg, f = (int(v) for v in d.ask('c %d %d %d %d' % (n, m, B1, mult * B1)))
            if sp is None or sq is None:
                degen += 1
                continue
            if sp[0] != sq[0]:
                exp = (1, p if sp[0] else q)
            else:
                w = uo.stage2_outcome(sp[1], sq[1])
                exp = (2, p) if w == 'p' else (2, q) if w == 'q' else (0, 1)
            ok = (stg, f) == exp
            if not ok and stg == 1 and f in (p, q):
                X, Z = (int(v) for v in d.ask('s %d %d %d' % (n, m, B1)))
                if X % f == 0 and Z % f == 0:
                    exc += 1
                    continue
            bad += not ok
            tot += 1
            stats[(B1, exp[0])] = stats.get((B1, exp[0]), 0) + 1
            if not ok and bad <= 8:
                print('mismatch B1 %d p %d q %d m %d: got stage %d f %d, expected %s'
                      % (B1, p, q, m, stg, f, exp))
    print('v2 curves: %d of %d mismatch (%d degenerate builds skipped, %d extra stage-1 finds '
          'from Edwards exceptions)' % (bad, tot, degen, exc))
    for B1, _, _ in em.LEVELS2:
        print('  B1 %3d: expected stage 1 %4d, stage 2 %4d, none %4d' % (
            B1, stats.get((B1, 1), 0), stats.get((B1, 2), 0), stats.get((B1, 0), 0)))
    return bad


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == 'curves':
        sys.exit(check(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 200,
                       int(sys.argv[4]) if len(sys.argv) > 4 else 1) != 0)
    print(__doc__)
