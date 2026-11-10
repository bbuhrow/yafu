#!/usr/bin/env python3
"""
Stage-1 scalars and chains for microecm / tinyecm (B1 = 27 .. 205).

Phase 0 of arith/micro_tiny_edwards_plan.md.

* Replays the existing x-only stage 1 exactly (uecm_stage1 / ecm_stage1:
  doublings while q < 4 B1, then the PRAC calls and the step tables of
  uprac70 / uprac85, parsed from microecm.c), on integer multiples.  That
  gives the scalar s of every level and the PRAC operation counts.
* Searches Edwards chains for the same s: a minimal-cost signed-digit
  representation for each odd digit set {1, 3, ..., m} (exact DP over the
  reachable values), plus the table cost.
* Prints the costs in multiplications (S = M, microecm) and with
  S = 0.8 M (tinyecm 2-word), and can emit the chains as a C header.

usage:
  tiny_chains.py                 report
  tiny_chains.py emit FILE.h     also write the digit tables
"""
import os, re, sys
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'factor', 'gmp-ecm')
LEVELS = (27, 47, 70, 85, 125, 165, 205)

# ---------------------------------------------------------------- PRAC replay

def parse_steps(fname, func):
    txt = open(os.path.join(SRC, fname)).read()
    i = txt.index(func)
    j = txt.index('{', txt.index('steps[', i))
    k = txt.index('}', j)
    return [int(v) for v in re.findall(r'\d+', txt[j:k])]


class Chain:
    """integer-multiple model of the x-only code; counts adds and dups"""
    def __init__(self):
        self.adds = 0
        self.dups = 0

    def add(self, a, b, c):
        assert abs(a - b) == c, (a, b, c)
        self.adds += 1
        return a + b

    def dup(self, a):
        self.dups += 1
        return 2 * a


def uprac(ch, P, c, v):
    shift = (c & -c).bit_length() - 1
    c >>= shift
    d = c
    r = int(float(d) * v + 0.5)
    d = c - r
    e = 2 * r - c
    p1 = p2 = p3 = P
    p1 = ch.dup(p1)
    while d != e:
        if d < e:
            d, e = e, d
            p1, p2 = p2, p1
        if d - e <= e // 4 and (d + e) % 3 == 0:
            d = (2 * d - e) // 3
            e = (e - d) // 2
            p4 = ch.add(p1, p2, p3)
            p5 = ch.add(p4, p1, p2)
            p2 = ch.add(p2, p4, p1)
            p1, p5 = p5, p1
        elif d - e <= e // 4 and (d - e) % 6 == 0:
            d = (d - e) // 2
            p2 = ch.add(p1, p2, p3)
            p1 = ch.dup(p1)
        elif (d + 3) // 4 <= e:
            d -= e
            p4 = ch.add(p2, p1, p3)
            p2, p4, p3 = p4, p3, p2
        elif (d + e) % 2 == 0:
            d = (d - e) // 2
            p2 = ch.add(p2, p1, p3)
            p1 = ch.dup(p1)
        elif d % 2 == 0:
            d //= 2
            p3 = ch.add(p3, p1, p2)
            p1 = ch.dup(p1)
        elif d % 3 == 0:
            d = d // 3 - e
            p4 = ch.dup(p1)
            p5 = ch.add(p1, p2, p3)
            p1 = ch.add(p4, p1, p1)
            p4 = ch.add(p4, p5, p3)
            p3, p2, p4 = p2, p4, p3
        elif (d + e) % 3 == 0:
            d = (d - 2 * e) // 3
            p4 = ch.add(p1, p2, p3)
            p2 = ch.add(p4, p1, p2)
            p4 = ch.dup(p1)
            p1 = ch.add(p1, p4, p1)
        elif (d - e) % 3 == 0:
            d = (d - e) // 3
            p4 = ch.add(p1, p2, p3)
            p3 = ch.add(p3, p1, p2)
            p2, p4 = p4, p2
            p4 = ch.dup(p1)
            p1 = ch.add(p1, p4, p1)
        else:
            e //= 2
            p3 = ch.add(p3, p2, p1)
            p2 = ch.dup(p2)
    P = ch.add(p1, p2, p3)
    for _ in range(shift):
        P = ch.dup(P)
    return P


def uprac_steps(ch, P, steps):
    p1 = p2 = p3 = P
    for st in steps:
        if st == 0:
            p1 = p2 = p3 = P
            p1 = ch.dup(p1)
        elif st == 3:
            p4 = ch.add(p2, p1, p3)
            p1, p4, p3, p2 = p4, p3, p2, p1
        elif st == 4:
            p4 = ch.add(p2, p1, p3)
            p2, p4, p3 = p4, p3, p2
        elif st == 5:
            p2 = ch.add(p2, p1, p3)
            p1 = ch.dup(p1)
        elif st == 6:
            P = ch.add(p1, p2, p3)
    return P


G = 0.618033988749894903
CALLS = {
    27: [(3, G), (3, G), (3, G), (5, G), (5, G), (7, G), (11, 0.580178728295464130),
         (13, G), (17, G), (19, G), (23, 0.522786351415446049)],
    47: [(3, G)] * 4 + [(5, G), (5, G), (7, G), (11, 0.580178728295464130), (13, G),
         (17, G), (19, G), (23, 0.522786351415446049), (29, 0.548409048446403258),
         (31, G), (37, 0.580178728295464130), (41, 0.548409048446403258), (43, G)],
}
_TAIL125 = [(5, G), (11, 0.580178728295464130), (89, G), (97, 0.723606797749978936),
            (101, 0.556250337855490828), (107, 0.580178728295464130),
            (109, 0.548409048446403258), (113, G)]
_TAIL165 = [(7747, 0.552188778811121), (131, G), (14111, 0.632839806088706),
            (20989, 0.620181980807415), (157, 0.640157392785047019),
            (163, 0.551390822543526449)]
CALLS[85] = ['p85', (61, 0.522786351415446049)]
CALLS[125] = ['p85'] + _TAIL125 + [(61, 0.522786351415446049), (103, 0.632839806088706269)]
CALLS[165] = ['p85'] + _TAIL125 + _TAIL165 + [(149, 0.580178728295464130)]
CALLS[205] = ['p85'] + _TAIL125 + _TAIL165 + [
    (13, G), (167, 0.580178728295464130), (173, 0.612429949509495031), (179, G),
    (181, 0.551390822543526449), (191, G), (193, G), (29353, 0.580178728295464),
    (199, 0.551390822543526449)]
CALLS[70] = ['p70']


def prac_level(B1, steps70, steps85):
    ch = Chain()
    P = 1
    q = 2
    while q < 4 * B1:
        P = ch.dup(P)
        q *= 2
    for c in CALLS[B1]:
        if c == 'p70':
            P = uprac_steps(ch, P, steps70)
        elif c == 'p85':
            P = uprac_steps(ch, P, steps85)
        else:
            P = uprac(ch, P, c[0], c[1])
    return P, ch.adds, ch.dups


def lcm_upto(B):
    from math import gcd
    l = 1
    for k in range(2, B + 1):
        l = l * k // gcd(l, k)
    return l

# ---------------------------------------------------------- Edwards chains

def edw_costs(sq):
    """M-equivalents of the extended a = -1 operations"""
    # mdup: Montgomery x-only doubling (3M + 2S); the power of 2 in s is
    # done after the hand-off to Montgomery (X:Z), where it is cheaper.
    return dict(dbl=3 + 4 * sq, dblT=4 + 4 * sq, add=7.0, cache=4.0, mdup=3 + 2 * sq)


def table_cost(m, C):
    """odd multiples 1..m in cached form (and extended form for the top)"""
    if m == 1:
        return C['cache']
    k = (m + 1) // 2
    # c(P), 2P with T, c(2P), then (k-1) x [add with T, cache]
    return C['cache'] + C['dblT'] + C['cache'] + (k - 1) * (C['add'] + 1 + C['cache'])


def best_repr(s, m, C):
    """min-cost chain for s from digits +-{1..m odd}: returns (cost, ops)"""
    D = list(range(1, m + 1, 2))
    addc = C['add'] + (C['dblT'] - C['dbl'])     # an add costs T in the dbl before it
    dblc = C['dbl']
    import sys as _s
    _s.setrecursionlimit(10000)

    @lru_cache(maxsize=None)
    def f(v):
        if v in Dset:
            return (0.0, ())
        if v % 2 == 0:
            c, ops = f(v // 2)
            return (c + dblc, ops + ('D',))
        best = None
        for t in D:
            for sgn in (1, -1):
                u = v - sgn * t
                if u <= 0:
                    continue
                c, ops = f(u)
                c += addc
                if best is None or c < best[0]:
                    best = (c, ops + (sgn * t,))
        return best

    Dset = frozenset(D)
    c, ops = f(s)
    f.cache_clear()
    return c, ops


def chain_eval(ops, top):
    """check a chain: start at top, D = double, int = add digit"""
    v = top
    for o in ops:
        v = 2 * v if o == 'D' else v + o
    return v


def edw_level(s, sq, mmax=31):
    C = edw_costs(sq)
    tz = (s & -s).bit_length() - 1
    so = s >> tz
    res = []
    for m in range(1, mmax + 1, 2):
        c, ops = best_repr(so, m, C)
        tc = table_cost(m, C)
        total = c + tc + tz * C['mdup']
        nd = sum(1 for o in ops if o == 'D') + tz
        na = sum(1 for o in ops if o != 'D')
        res.append((total, m, nd, na, tc, ops, tz))
    res.sort(key=lambda r: r[0])
    return res


def top_of(s, ops, tz):
    """the starting digit of a chain (ops are applied after it)"""
    v = s >> tz
    for o in reversed(ops):
        v = v // 2 if o == 'D' else v - o
    return v


def main():
    steps70 = parse_steps('microecm.c', 'static void uecm_uprac70(')
    steps85 = parse_steps('microecm.c', 'static void uecm_uprac85(')
    t70 = parse_steps('tinyecm.c', 'void prac70(')
    t85 = parse_steps('tinyecm.c', 'void prac85(')
    print('tinyecm step tables identical to microecm: 70 %s, 85 %s' %
          (t70 == steps70, t85 == steps85))
    emit = len(sys.argv) > 2 and sys.argv[1] == 'emit'
    out = []
    print('\nlevel  bits | PRAC adds dups   M(S=M) M(S=.8) | '
          'Edw m  dbls adds table   M(S=M) ratio | m(S=.8)   M(S=.8) ratio')
    for B1 in LEVELS:
        s, adds, dups = prac_level(B1, steps70, steps85)
        # sanity: s is B1-powersmooth-ish: every prime < B1 divides it
        lc = lcm_upto(B1 - 1)
        missing = [q for q in range(2, B1) if all(q % d for d in range(2, q)) and s % q]
        bits = s.bit_length()
        pm = adds * 6 + dups * 5
        pm8 = adds * (4 + 2 * .8) + dups * (3 + 2 * .8)
        r1 = edw_level(s, 1.0)
        r8 = edw_level(s, 0.8)
        b1, b8 = r1[0], r8[0]
        print('%5d %5d | %9d %4d %8.0f %7.0f | %6d %5d %4d %5.0f %8.0f %5.3f | %5d %9.0f %5.3f%s'
              % (B1, bits, adds, dups, pm, pm8, b1[1], b1[2], b1[3], b1[4], b1[0], pm / b1[0],
                 b8[1], b8[0], pm8 / b8[0], ('  missing primes %s' % missing) if missing else ''))
        assert chain_eval(b1[5], top_of(s, b1[5], b1[6])) << b1[6] == s
        out.append((B1, s, b1))
        if '-v' in sys.argv:
            for r in r1[:6]:
                print('        m=%2d  M=%.0f dbls=%d adds=%d table=%.0f' % (r[1], r[0], r[2], r[3], r[4]))
    if emit:
        write_header(sys.argv[2], out)


def write_header(fname, out):
    L = []
    L.append('/* generated by arith/edwards/tiny_chains.py: Edwards stage-1 chains for')
    L.append('   the microecm / tinyecm levels.  Each chain: table of odd multiples')
    L.append('   1..m, start at digit "top", then ops: 0 = double, +-k = add +-kP.')
    L.append('   Then "tz" Montgomery x-only doublings after the hand-off.')
    L.append('   s matches uecm_stage1 exactly. */')
    L.append('#ifndef EDW_TINY_CHAINS_H')
    L.append('#define EDW_TINY_CHAINS_H')
    L.append('#include <stdint.h>')
    for B1, s, b in out:
        total, m, nd, na, tc, ops, tz = b
        top = top_of(s, ops, tz)
        enc = [0 if o == 'D' else o for o in ops]
        L.append('static const int8_t edw_chain_%d[%d] = {' % (B1, len(enc)))
        for i in range(0, len(enc), 24):
            L.append('    ' + ','.join(str(v) for v in enc[i:i + 24]) + ',')
        L.append('};')
        L.append('#define EDW_CHAIN_%d_LEN %d' % (B1, len(enc)))
        L.append('#define EDW_CHAIN_%d_M %d' % (B1, m))
        L.append('#define EDW_CHAIN_%d_TOP %d' % (B1, top))
        L.append('#define EDW_CHAIN_%d_TZ %d' % (B1, tz))
        L.append('#define EDW_CHAIN_%d_ADDS %d' % (B1, na))
    L.append('#endif')
    open(fname, 'w', newline='\n').write('\n'.join(L) + '\n')
    print('wrote', fname)


if __name__ == '__main__':
    main()
