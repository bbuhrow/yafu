#!/usr/bin/env python3
"""
Small-height Edwards curves for microecm / tinyecm, and the yield comparison
at their sizes (Phase 0 of arith/micro_tiny_edwards_plan.md).

Curve families (a x^2 + y^2 = 1 + d x^2 y^2, d = dn/dd):
  iv    a = -1, d = -e^4, e = (g - 1/g)/2, g = gn/gd   (GKL family (iv):
        e^2 + 1 is a square; Z/2 x Z/4 torsion, 16 | #E)
  z12   a = +1, EECM Thm 4.9: d = (u^2+1)^3 (u^2-4u+1) / ((u-1)^6 (u+1)^2)
  z2x8  a = +1, EECM Thm 4.5: x8 = (u^2+2u+2)/(u^2-2), d = (2 x8^2 - 1)/x8^4
For each parameter a small non-torsion point is searched (x = xn/xd, then
y^2 = (1 - a x^2)/(1 - d x^2) must be a rational square).  Kept: curves
whose projective start point X = xn yd, Y = yn xd, Z = xd yd, T = xn yn and
dn, dd all fit in 63 bits.  The table is sorted by height.

usage:
  tiny_curves.py search FAM [maxparam] [maxpt] [maxbits]   -> curves_FAM.json
  tiny_curves.py yield pbits B1 nprimes [seed]    paired yield at one level
  tiny_curves.py yieldall nprimes [seed]          all levels
"""
import json
import math
import os
import random
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'vecpoly', 'ecm_tests'))
import bs_oracle as bo                      # noqa: E402
import tiny_chains as tc                    # noqa: E402
import edw_ref                              # noqa: E402

LIM = 1 << 63

# ----------------------------------------------------------------- curves

def ed_add_q(P, Q, a, d):
    (x1, y1), (x2, y2) = P, Q
    t = d * x1 * x2 * y1 * y2
    if t == 1 or t == -1:
        return None
    return ((x1 * y2 + y1 * x2) / (1 + t), (y1 * y2 - a * x1 * x2) / (1 - t))


def is_torsion(P, a, d, kmax=16):
    if P[0] == 0 or P[1] == 0:
        return True
    R = P
    for _ in range(2, kmax + 1):
        R = ed_add_q(R, P, a, d)
        if R is None or R == (0, 1):
            return True
        # points of order dividing 2k: small x or y signals torsion too
    return False


def params(fam, maxparam):
    out = []
    seen = set()
    for den in range(1, maxparam + 1):
        for num in range(-maxparam, maxparam + 1):
            if math.gcd(num, den) != 1:
                continue
            u = F(num, den)
            try:
                if fam == 'iv':
                    if u <= 0 or u == 1:
                        continue
                    e = (u - 1 / u) / 2
                    a, d = -1, -e ** 4
                elif fam == 'z12':
                    if u in (0, 1, -1):
                        continue
                    a = 1
                    d = (u * u + 1) ** 3 * (u * u - 4 * u + 1) / ((u - 1) ** 6 * (u + 1) ** 2)
                elif fam == 'z2x8':
                    a = 1
                    if u * u == 2:
                        continue
                    x8 = (u * u + 2 * u + 2) / (u * u - 2)
                    if x8 == 0:
                        continue
                    d = (2 * x8 * x8 - 1) / x8 ** 4
                else:
                    raise ValueError(fam)
            except ZeroDivisionError:
                continue
            if d in (0, 1, -1) or d == a or d in seen:
                continue
            seen.add(d)
            out.append((u, a, d))
    return out


def find_point(a, d, maxpt):
    """smallest-height non-torsion point with x = xn/xd, xn, xd <= maxpt"""
    dn, dd = d.numerator, d.denominator
    best = None
    for s in range(2, 2 * maxpt + 1):           # by xn + xd
        for xd in range(1, s):
            xn = s - xd
            if xn > maxpt or xd > maxpt or math.gcd(xn, xd) != 1:
                continue
            # y^2 = (xd^2 - a xn^2) / (xd^2 - d xn^2)
            num = xd * xd - a * xn * xn
            den = dd * xd * xd - dn * xn * xn          # times dd
            if den == 0 or num == 0:
                continue
            w = num * dd * den
            if w <= 0:
                continue
            r = math.isqrt(w)
            if r * r != w:
                continue
            y = F(r, den)                               # y^2 = num dd / den
            if y < 0:
                y = -y
            x = F(xn, xd)
            assert a * x * x + y * y == 1 + d * x * x * y * y
            if is_torsion((x, y), a, d):
                continue
            X = x.numerator * y.denominator
            Y = y.numerator * x.denominator
            Z = x.denominator * y.denominator
            T = x.numerator * y.numerator
            h = max(X, abs(Y), Z, abs(T), abs(dn), dd)
            if best is None or h < best[0]:
                best = (h, x, y)
        if best is not None and s > 40:
            break
    return best


def _search_one(job):
    fam, u, a, d, maxpt = job
    b = find_point(a, d, maxpt)
    if b is None:
        return None
    h, x, y = b
    return dict(fam=fam, param=str(u), a=a, dn=d.numerator, dd=d.denominator,
                x=str(x), y=str(y), height_bits=h.bit_length())


def search(fam, maxparam=24, maxpt=200, maxbits=63):
    import multiprocessing as mp
    jobs = [(fam, u, a, d, maxpt) for (u, a, d) in params(fam, maxparam)]
    with mp.Pool() as pool:
        res = [r for r in pool.map(_search_one, jobs, chunksize=4) if r]
    res = [r for r in res if r['height_bits'] <= maxbits]
    res.sort(key=lambda r: r['height_bits'])
    fn = os.path.join(HERE, 'curves_%s.json' % fam)
    json.dump(res, open(fn, 'w'), indent=1)
    print('%s: %d parameters, %d curves with a point of height <= %d bits -> %s'
          % (fam, len(jobs), len(res), maxbits, fn))
    hist = {}
    for r in res:
        hist[r['height_bits']] = hist.get(r['height_bits'], 0) + 1
    print('height histogram:', dict(sorted(hist.items())))
    return res


def load(fam):
    return json.load(open(os.path.join(HERE, 'curves_%s.json' % fam)))

# ------------------------------------------------------------------ yield

def inv(v, p):
    return pow(v % p, -1, p)


def table_curve_point(c, p):
    """Montgomery image (A, B, (u, v)) of a table curve mod p, or None if
    degenerate mod p (a build "find")"""
    a = c['a']
    dn, dd = c['dn'], c['dd']
    x, y = F(c['x']), F(c['y'])
    try:
        d = dn * inv(dd, p) % p
        xm = x.numerator * inv(x.denominator, p) % p
        ym = y.numerator * inv(y.denominator, p) % p
        if (a - d) % p == 0 or (1 - ym) % p == 0 or xm == 0:
            return None
        A = 2 * (a + d) * inv(a - d, p) % p
        B = 4 * inv(a - d, p) % p
        uu = (1 + ym) * inv(1 - ym, p) % p
        vv = uu * inv(xm, p) % p
        return A, B, (uu, vv)
    except ValueError:
        return None


def stage2_pairs(B1):
    """(A, b) pairs of microecm's stage 2 at this level: a prime r is found
    when r | A - b or r | A + b"""
    src = open(os.path.join(tc.SRC, 'microecm.c')).read()

    def arr(name):
        i = src.index('uint8_t %s[' % name)
        j = src.index('{', i)
        k = src.index('}', j)
        return [int(v) for v in src[j + 1:k].replace('\n', ' ').split(',') if v.strip()]

    pairs = []
    if B1 < 70:
        st = [59, 53, 49, 47, 43, 37, 31, 29, 23, 19, 17, 11, 7, 1, 13, 41] if B1 == 27 else \
             [43, 37, 31, 29, 23, 19, 17, 11, 7, 1, 13, 41, 47, 49, 59]
        pairs += [(90, b) for b in st]
        Pa = 150
    elif B1 == 70:
        st = [49, 47, 41, 37, 31, 23, 19, 17, 13, 11, 7, 29, 43, 53, 59]
        pairs += [(120, b) for b in st]
        Pa = 150
    elif B1 == 165:
        Pa = 180
    elif B1 == 205:
        Pa = 210
    else:
        Pa = 150
    name = {27: 'b1_27', 47: 'b1_47', 70: 'b1_70', 85: 'b1_85',
            125: 'b1_125_x50', 165: 'b1_165_x50', 205: 'b1_205_x50'}[B1]
    ba = arr(name)
    i = 0
    while i < len(ba):
        if ba[i] == 0:
            Pa += 120
            i += 1
        pairs.append((Pa, ba[i]))
        i += 1
    return pairs


_LEVEL = {}


def level_data(B1):
    if B1 not in _LEVEL:
        s70 = tc.parse_steps('microecm.c', 'static void uecm_uprac70(')
        s85 = tc.parse_steps('microecm.c', 'static void uecm_uprac85(')
        s, _, _ = tc.prac_level(B1, s70, s85)
        pr = stage2_pairs(B1)
        vals = set()
        for A, b in pr:
            vals.add(A - b)
            vals.add(A + b)
        _LEVEL[B1] = (s, sorted(vals))
    return _LEVEL[B1]


def outcome(A, g, P, p, s, vals):
    r0 = bo.point_order(P, A, g, p)
    r = r0 // math.gcd(r0, s)
    if r == 1:
        return 1, r0
    if any(v % r == 0 for v in vals):
        return 2, r0
    return 0, r0


def _yield_one(args):
    p, B1, picks = args
    s, vals = level_data(B1)
    out = {}
    for fam, c in picks.items():
        if fam == 'suyama':
            x0, A, g = bo.suyama(c, p)
            if g == 0:
                out[fam] = (1, 0)
                continue
            st, r0 = outcome(A, g, (x0, 1), p, s, vals)
        elif fam == 'ivr':                      # GKL family (iv), random m
            try:
                A, g, P = edw_ref.curve_point('iv', c, p)
            except (edw_ref.Find, ValueError, ZeroDivisionError):
                out[fam] = (1, 0)
                continue
            st, r0 = outcome(A, g, P, p, s, vals)
        else:
            cp = table_curve_point(c, p)
            if cp is None:
                out[fam] = (1, 0)
                continue
            A, B, P = cp
            st, r0 = outcome(A, B, P, p, s, vals)
        out[fam] = (st, tc_v2(r0))
    return out


def tc_v2(n):
    c = 0
    while n and n % 2 == 0:
        n //= 2
        c += 1
    return c


FAMS = ('iv', 'z12', 'z2x8')


def yield_cmp(pbits, B1, nprimes, seed=1, tables=None, K=64):
    import multiprocessing as mp
    rng = random.Random(seed * 1000003 + pbits * 1009 + B1)
    tables = tables or {f: load(f)[:K] for f in FAMS}
    jobs = []
    for _ in range(nprimes):
        p = rng.randrange(1 << (pbits - 1), 1 << pbits) | 1
        while not bo.is_prime(p):
            p += 2
        picks = {'suyama': rng.randrange(7, 1 << 32), 'ivr': rng.randrange(1, 1 << 31)}
        for f in FAMS:
            picks[f] = rng.choice(tables[f])
        jobs.append((p, B1, picks))
    with mp.Pool() as pool:
        res = pool.map(_yield_one, jobs, chunksize=64)
    row = {}
    for f in ('suyama', 'ivr') + FAMS:
        row[f] = (sum(1 for r in res if r[f][0] == 1), sum(1 for r in res if r[f][0] >= 1),
                  sum(r[f][1] for r in res) / len(res))
    base = row['suyama'][1]
    line = '%3d %4d %6d | ' % (pbits, B1, nprimes)
    for f in ('suyama', 'ivr') + FAMS:
        s1, s12, v2 = row[f]
        line += '%-6s %5d/%5d %.3f v2 %.2f | ' % (f, s1, s12, s12 / base, v2)
    # paired std error for each family vs suyama
    errs = []
    for f in ('ivr',) + FAMS:
        dif = [int(r[f][0] >= 1) - int(r['suyama'][0] >= 1) for r in res]
        mu = sum(dif) / len(dif)
        sd = math.sqrt(sum((x - mu) ** 2 for x in dif) / (len(dif) - 1)) * math.sqrt(len(dif))
        errs.append('%s %+d+-%.0f' % (f, sum(dif), sd))
    print(line + ' '.join(errs), flush=True)
    return row


# dispatch levels: (factor bits, B1), as microecm / tinyecm pick them
LEVELS = [(20, 27), (22, 47), (24, 70), (26, 85), (29, 125), (31, 165), (32, 205)]


def main():
    cmd = sys.argv[1]
    if cmd == 'search':
        fam = sys.argv[2]
        mp_ = int(sys.argv[3]) if len(sys.argv) > 3 else 24
        mx = int(sys.argv[4]) if len(sys.argv) > 4 else 200
        mb = int(sys.argv[5]) if len(sys.argv) > 5 else 63
        search(fam, mp_, mx, mb)
    elif cmd == 'yield':
        yield_cmp(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]),
                  int(sys.argv[5]) if len(sys.argv) > 5 else 1)
    elif cmd == 'yieldall':
        n = int(sys.argv[2])
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        for pb, B1 in LEVELS:
            yield_cmp(pb, B1, n, seed)
    else:
        print(__doc__)


if __name__ == '__main__':
    main()
