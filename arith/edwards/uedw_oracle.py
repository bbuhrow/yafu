#!/usr/bin/env python3
"""
Exact oracle check of factor/gmp-ecm/microecm_edw.c (Phase 1 of
arith/micro_tiny_edwards_plan.md), through the uedw_test driver.

  uedw_oracle.py params BIN [ncases]
      curve parameters for m = 1..40 mod random primes vs edw_ref.curve_iv
  uedw_oracle.py curves BIN [nper] [seed]
      one curve per (n = p q, m) at every level: the stage of the find vs
      the prediction from exact ord P mod p and mod q, microecm's s and its
      stage-2 pairs
"""
import math
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'vecpoly', 'ecm_tests'))
import bs_oracle as bo          # noqa: E402
import edw_ref as er            # noqa: E402
import tiny_curves as tcv       # noqa: E402


class Driver:
    def __init__(self, path):
        self.p = subprocess.Popen([path], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  universal_newlines=True, bufsize=1)

    def ask(self, line):
        self.p.stdin.write(line + '\n')
        self.p.stdin.flush()
        return self.p.stdout.readline().split()


def rand_prime(rng, bits):
    while True:
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if bo.is_prime(p):
            return p


TABLE = 0x80000000
_TAB = None


def table_curves():
    """the constants of factor/gmp-ecm/edw_tiny_curves.h: X, Y, Z, T, dn, dd"""
    global _TAB
    if _TAB is None:
        import re
        src = open(os.path.join(tcv.tc.SRC, 'edw_tiny_curves.h')).read()
        _TAB = [[int(v) for v in re.findall(r'-?\d+(?=LL)', row)]
                for row in re.findall(r'\{([^{}]*LL[^{}]*)\}', src)]
    return _TAB


def curve_ref(cid, p):
    """(d, x, y) mod p of curve id cid (GKL m, or TABLE | i); raises er.Find"""
    if cid & TABLE:
        X, Y, Z, T, dn, dd = table_curves()[(cid & ~TABLE) % len(table_curves())]
        zi = er.inv(Z, p)
        return dn * er.inv(dd, p) % p, X * zi % p, Y * zi % p
    return er.curve_iv(cid, p)


def check_params(binp, ncases=20):
    d = Driver(binp)
    rng = random.Random(5)
    bad = tot = 0
    for _ in range(ncases):
        p = rand_prime(rng, rng.choice([40, 52, 61, 63, 64]))
        ntab = len(table_curves())
        for m in list(range(1, 41)) + [TABLE | i for i in range(ntab)]:
            X, Y, Z, T, dn, dd = (int(v) for v in d.ask('p %d %d' % (p, m)))
            try:
                dref, x1, y1 = curve_ref(m, p)
            except er.Find:
                continue
            zi = pow(Z, -1, p)
            ok = (X * zi % p == x1 and Y * zi % p == y1 and T * zi % p == x1 * y1 % p
                  and dn * pow(dd, -1, p) % p == dref)
            bad += not ok
            tot += 1
            if not ok and bad < 5:
                print('params mismatch p %d m %d' % (p, m))
    print('params: %d of %d mismatch' % (bad, tot))
    return bad


MAP = [0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
       0, 0, 0, 8, 0, 0, 0, 0, 0, 9, 0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
       0, 12, 0, 13, 0, 0, 0, 14, 0, 15, 0, 0, 0, 16, 0, 0, 0, 0, 0, 17]
_BARRAY = {}


def barray(B1):
    if B1 not in _BARRAY:
        src = open(os.path.join(tcv.tc.SRC, 'microecm.c')).read()
        name = {27: 'b1_27', 47: 'b1_47', 70: 'b1_70', 85: 'b1_85', 125: 'b1_125_x50',
                165: 'b1_165_x50', 205: 'b1_205_x50'}[B1]
        i0 = src.index('uint8_t %s[' % name)
        body = src[src.index('{', i0) + 1:src.index('}', i0)]
        _BARRAY[B1] = [int(v) for v in body.split(',') if v.strip()]
    return _BARRAY[B1]


def stage2_zero_flags(P0, A, p, B1):
    """microecm's x-only stage 2, exactly, mod p from the stage-1 point P0 =
    (X, Z): for every cross term, whether it is 0 mod p, as two lists (the
    B1 <= 70 pre-steps, and the main pairing loop)."""
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
    pt90, pt60 = Pa, Pd
    Pd = uadd(Pa, Pad, Pd)
    Pa = uadd(Pd, Pad, Pa)

    def cross(R, b):
        Rprod = R[0] * R[1] % p
        tt1 = (R[0] - Pb[MAP[b]][0]) * (R[1] + Pb[MAP[b]][1]) % p
        return (tt1 + Pbprod[MAP[b]] - Rprod) % p

    pre, main = [], []
    if B1 < 70:
        steps = ([59, 53, 49, 47, 43, 37, 31, 29, 23, 19, 17, 11, 7, 1, 13, 41] if B1 == 27
                 else [43, 37, 31, 29, 23, 19, 17, 11, 7, 1, 13, 41, 47, 49, 59])
        for b in steps:
            pre.append(cross(pt90, b) == 0)
    elif B1 == 70:
        for b in [49, 47, 41, 37, 31, 23, 19, 17, 13, 11, 7, 29, 43, 53, 59]:
            pre.append(cross(Pd, b) == 0)
    elif B1 == 165:
        Pa = dup(pt90)
        Pad = pt60
    elif B1 == 205:
        Pa = uadd(Pd, pt90, Pad)
        Pad = pt90
    ba = barray(B1)
    i = 0
    while i < len(ba):
        if ba[i] == 0:
            point = Pa
            Pa = uadd(Pa, Pd, Pad)
            Pad = point
            i += 1
        main.append(cross(Pa, ba[i]) == 0)
        i += 1
    return pre, main


def stage2_outcome(fp, fq):
    """replay the C accumulator loops (each breaks when a term would make acc
    0 mod n = p q, keeping acc) on the zero flags of p and q; returns which
    prime acc is 0 mod: 'p', 'q' or None"""
    zp = zq = False
    for k in (0, 1):
        for a, b in zip(fp[k], fq[k]):
            np_, nq = zp or a, zq or b
            if np_ and nq:
                break
            zp, zq = np_, nq
    if zp != zq:
        return 'p' if zp else 'q'
    return None


def predict(p, m, B1, s):
    """(stage-1 found, stage-2 zero flags) mod p, or None if the curve build
    is degenerate mod p"""
    try:
        d, x, y = curve_ref(m, p)
        if (1 + d) % p == 0:
            raise er.Find()
        A, g, u, v = er.edw_to_mont(d, x, y, p)
        P = (u, v)
    except (er.Find, ValueError, ZeroDivisionError):
        return None
    R0, _ = bo.x_ladder(s, P[0], A, p)
    return (R0[1] % p == 0, stage2_zero_flags(R0, A, p, B1))


LEVELS = [(20, 27), (22, 47), (24, 70), (26, 85), (29, 125), (31, 165), (32, 205)]


def check_curves(binp, nper=300, seed=1):
    d = Driver(binp)
    rng = random.Random(seed)
    tot = bad = degen = exc = 0
    stats = {}
    for pb, B1 in LEVELS:
        s, _ = tcv.level_data(B1)
        for _ in range(nper):
            p = rand_prime(rng, pb)
            q = rand_prime(rng, pb + rng.randrange(0, 4))
            if p == q:
                continue
            n = p * q
            if n >= 1 << 64:
                continue
            m = rng.randrange(1, 400)
            if rng.random() < 0.3:
                m = TABLE | rng.randrange(len(table_curves()))
            sp, sq = predict(p, m, B1, s), predict(q, m, B1, s)
            st, f = (int(v) for v in d.ask('c %d %d %d %d' % (n, m, B1, 25 * B1)))
            if sp is None or sq is None:
                degen += 1          # curve build degenerate mod p or q
                continue
            if sp[0] != sq[0]:
                exp = (1, p if sp[0] else q)
            else:
                w = stage2_outcome(sp[1], sq[1])
                exp = (2, p) if w == 'p' else (2, q) if w == 'q' else (0, 1)
            ok = (st, f) == exp
            if not ok and st == 1 and f in (p, q):
                # an exceptional case of the Edwards formulas mod f: the
                # stage-1 point is (0 : 0) mod f, a correct "free" find
                X, Z = (int(v) for v in d.ask('s %d %d %d' % (n, m, B1)))
                if X % f == 0 and Z % f == 0:
                    exc += 1
                    continue
            bad += not ok
            tot += 1
            stats[(B1, exp[0])] = stats.get((B1, exp[0]), 0) + 1
            if not ok and bad <= 8:
                print('mismatch B1 %d p %d q %d m %d: got stage %d f %d, expected %s'
                      % (B1, p, q, m, st, f, exp))
    print('curves: %d of %d mismatch (%d degenerate builds skipped, %d extra stage-1 finds '
          'from Edwards exceptions)' % (bad, tot, degen, exc))
    for B1 in [b for _, b in LEVELS]:
        print('  B1 %3d: expected stage 1 %4d, stage 2 %4d, none %4d' % (
            B1, stats.get((B1, 1), 0), stats.get((B1, 2), 0), stats.get((B1, 0), 0)))
    return bad


if __name__ == '__main__':
    cmd, binp = sys.argv[1], sys.argv[2]
    if cmd == 'params':
        sys.exit(check_params(binp, int(sys.argv[3]) if len(sys.argv) > 3 else 20) != 0)
    elif cmd == 'curves':
        sys.exit(check_curves(binp, int(sys.argv[3]) if len(sys.argv) > 3 else 300,
                              int(sys.argv[4]) if len(sys.argv) > 4 else 1) != 0)
    else:
        print(__doc__)
