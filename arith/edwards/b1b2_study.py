#!/usr/bin/env python3
"""
B1 / B2 study for microecm_edw (arith/micro_tiny_edwards_plan.md, F5 / Part I).

For inputs n with factors of N bits, which stage-1 bound B1 (and scalar s)
and which stage-2 bound B2 minimize the expected time to a factor?

  b1b2_study.py gen N0 N1 NP K [seed]
      exact point orders: for NP random primes of each size N0..N1, the
      order of the start point of curves 0..K-1 of microecm_edw's sequence
      (table curves 0..4, then GKL [m]Q, m = m0, m0 + 1, ..., m0 random in
      [3, 258] per prime).  Writes b1b2_orders_N.txt (one line per prime:
      p then K orders; 0 = degenerate build, a find).
  b1b2_study.py eval N0 N1 MODE [out.txt]
  b1b2_study.py decomp N0 N1 MODE [out.txt]
  b1b2_study.py plan N0 N1 MODE B1,B1,... [out.txt]
      MODE = one: n has one N-bit factor (the cofactor is never found);
      MODE = two: n = p q with p, q both N bits (balanced; consecutive rows
      of the order file pair up).  eval: the best candidates per N and the
      flatness; decomp: which change gives what; plan: one (s policy, B2
      multiplier) per level and a level per N, vs the legacy dispatch.

Success model per curve (exact, as uedw_oracle.py checks it): with
r = ord / gcd(ord, s), stage 1 finds p when r = 1; stage 2 when r divides
one of the values A - b, A + b of the stage-2 table (microecm's own tables
for the legacy levels, s2_tables.table for others).  So stage 2 also
catches extra powers of primes >= 7 and some composite r, never extra
powers of 2, 3 or 5 (A +- b is coprime to 30).

Cost model (ns per curve): calibrated with uedw_time.c, then scaled by CAL
to the realistic per-curve times of uedw_lvlbench.c (varying inputs and
curves cost ~1.33x the repeated-input microbenchmark, evenly across B1).
"""
import math
import os
import random
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'vecpoly', 'ecm_tests'))
import bs_oracle as bo          # noqa: E402
import edw_ref as er            # noqa: E402
import tiny_chains as tc        # noqa: E402
import tiny_curves as tcv       # noqa: E402
import s2_tables as s2          # noqa: E402
import uedw_oracle as uo        # noqa: E402

NTAB = 5


# ------------------------------------------------------------------ orders

def curve_mont(k, m0, p):
    """(A, B, P) of curve k of the sequence mod p; None if degenerate"""
    try:
        if k < NTAB:
            d, x, y = uo.curve_ref(uo.TABLE | k, p)
        else:
            d, x, y = er.curve_iv(m0 + k - NTAB, p)
        if (1 + d) % p == 0:
            return None
        A, B, u, v = er.edw_to_mont(d, x, y, p)
        return A, B, (u, v)
    except (er.Find, ValueError, ZeroDivisionError):
        return None


def _orders(args):
    p, m0, K = args
    out = []
    for k in range(K):
        c = curve_mont(k, m0, p)
        if c is None:
            out.append(0)
            continue
        A, B, P = c
        try:
            out.append(bo.point_order(P, A, B, p))
        except (AssertionError, ValueError, ZeroDivisionError):
            out.append(0)
    return p, out


def gen(N0, N1, NP, K, seed=1):
    rng = random.Random(seed)
    with Pool() as pool:
        for N in range(N0, N1 + 1):
            jobs = []
            seen = set()
            while len(jobs) < NP:
                p = uo.rand_prime(rng, N)
                if p in seen or p < 7:
                    continue
                seen.add(p)
                jobs.append((p, rng.randrange(3, 259), K))
            res = pool.map(_orders, jobs, chunksize=4)
            with open(os.path.join(HERE, 'b1b2_orders_%d.txt' % N), 'w') as f:
                for p, o in res:
                    f.write('%d %s\n' % (p, ' '.join(map(str, o))))
            print('N %d: %d primes x %d curves' % (N, NP, K), flush=True)


def load(N):
    rows = []
    with open(os.path.join(HERE, 'b1b2_orders_%d.txt' % N)) as f:
        for line in f:
            v = [int(t) for t in line.split()]
            rows.append((v[0], v[1:]))
    return rows


# ------------------------------------------------------------- scalars, costs

def primes_upto(n):
    return s2.primes_upto(n)


def scalar(B1, pol):
    """stage-1 scalar.  pol = (k2, kodd): q^e <= k B1 (k2 for q = 2,
    kodd for odd q), every prime q <= B1."""
    k2, ko = pol
    s = 1
    for q in primes_upto(B1):
        lim = (k2 if q == 2 else ko) * B1
        e = q
        while e * q <= lim:
            e *= q
        s *= e
    return s


_LEVEL_S = {}


def level_scalar(B1):
    """microecm's exact s at an existing level"""
    if B1 not in _LEVEL_S:
        s70 = tc.parse_steps('microecm.c', 'static void uecm_uprac70(')
        s85 = tc.parse_steps('microecm.c', 'static void uecm_uprac85(')
        _LEVEL_S[B1] = tc.prac_level(B1, s70, s85)[0]
    return _LEVEL_S[B1]


def get_scalar(B1, pol):
    return level_scalar(B1) if pol == 'lvl' else scalar(B1, pol)


_CHAIN = {}


def chain_m(s):
    """M-equivalents of the best Edwards chain for s (tiny_chains' DP)"""
    if s not in _CHAIN:
        _CHAIN[s] = tc.edw_level(s, 1.0, mmax=15)[0][0]
    return _CHAIN[s]


# cost model (ns), fitted to uedw_time.c on the i5-11600K (one core):
#   stage 1 = S1_FIX + S1_PER_M * chain M;  stage 2 = S2_FIX + S2_X * (crosses
#   + S2_G * giants);  per curve of the GKL sequence a build of ~46 M; per
#   input a context setup and, past the table curves, the [m0]Q start.
#   All times CAL x (see the module doc).
#   Then refit (b1b2_fit.py) to 45 uedw_lvlbench.c configurations (N 18-26,
#   legacy and v2 levels, realistic inputs): rms error 1.9%.  The refit
#   stage-2 cost per cross term is ~2x the microbenchmark's.
CAL = 1.33
S1_FIX, S1_PER_M = 0.0, 1.155
S2_FIX, S2_X, S2_G = 410.0, 3.299, 3.5
GKL_BUILD = 78.0
INPUT_FIX = 80.0
GKL_START = 200.0 * CAL

LEVELS = [27, 47, 70, 85, 125, 165, 205]


def legacy_mult(B1):
    return 50 if B1 >= 125 else 25


def t1_cost(B1, pol):
    return S1_FIX + S1_PER_M * chain_m(get_scalar(B1, pol))


class Stage2:
    """a stage-2 table: its cost, and which r it catches (r | A +- b)"""

    def __init__(self, B1, mult, pol):
        self.B2 = mult * B1
        self.cache = {}
        if mult == 0:
            self.cost = 0.0
            self.vals = []
        elif pol == 'lvl' and B1 in LEVELS and mult == legacy_mult(B1):
            pr = tcv.stage2_pairs(B1)               # microecm's own table
            self.vals = sorted(set(v for A, b in pr for v in (abs(A - b), A + b)))
            ng = len(set(A for A, _ in pr))
            self.cost = S2_FIX + S2_X * (len(pr) + S2_G * ng)
        else:
            tab = s2.table(B1, self.B2)
            self.vals = sorted(set(v for A, bs in tab for b in bs for v in (abs(A - b), A + b)))
            nc, ng = s2.cost(tab)
            ng += (tab[0][0] - 30) // 120 if tab else 0     # skipped windows from A = 30
            self.cost = S2_FIX + S2_X * (nc + S2_G * ng)
        self.vset = set(self.vals)
        self.mx = self.vals[-1] if self.vals else 0

    def catches(self, r):
        if r > self.mx:
            return False
        c = self.cache.get(r)
        if c is None:
            c = r in self.vset or any(v % r == 0 for v in self.vals)
            self.cache[r] = c
        return c


_S2 = {}


def stage2(B1, mult, pol):
    key = (B1, mult, pol == 'lvl' and B1 in LEVELS and mult == legacy_mult(B1))
    if key not in _S2:
        _S2[key] = Stage2(B1, mult, pol)
    return _S2[key]


# -------------------------------------------------------------------- eval

def curve_outcome(o, s, st2):
    """0 = none, 1 = stage 1, 2 = stage 2"""
    if o == 0:
        return 1
    r = o // math.gcd(o, s)
    if r == 1:
        return 1
    if st2.vals and st2.catches(r):
        return 2
    return 0


def eval_one(rows, s, st2, t1, mode, K=48):
    """time per input with a find (ns), fraction found within K curves,
    mean curves.  mode 'one': each row is an input; 'two': rows 2i, 2i+1
    are the two factors of one input (a curve finding both in stage 1
    finds nothing; in stage 2 the accumulator stops at the first)."""
    t2 = st2.cost
    tot = 0.0
    found = 0
    ncur = 0
    if mode == 'one':
        inputs = [(o,) for _, o in rows]
    else:
        inputs = [(rows[i][1], rows[i + 1][1]) for i in range(0, len(rows) - 1, 2)]
    for ords in inputs:
        t = INPUT_FIX
        for k in range(K):
            if k >= NTAB:
                t += GKL_BUILD + (GKL_START if k == NTAB else 0.0)
            t += t1
            outs = [curve_outcome(o[k], s, st2) for o in ords]
            n1 = outs.count(1)
            if n1 == 1:
                found += 1
                ncur += k + 1
                break
            if t2:
                t += t2
                if n1 == 0 and 2 in outs:
                    found += 1
                    ncur += k + 1
                    break
        else:
            ncur += K
        tot += t
    n = len(inputs)
    return tot / max(found, 1), found / n, ncur / n


B1_GRID = [15, 20, 25, 27, 30, 35, 40, 47, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 110, 125,
           140, 150, 165, 180, 205, 230, 260, 300, 350, 400]
MULTS = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 200]
POLS = [(1, 1), (4, 1), (4, 2), (16, 4), (2, 1), (8, 2)]


def candidates(B1s, pols=None):
    for B1 in B1s:
        for pol in (pols or POLS + (['lvl'] if B1 in LEVELS else [])):
            if pol == 'lvl' and B1 not in LEVELS:
                continue
            for m in MULTS:
                yield B1, pol, m


def ev(rows, B1, pol, m, mode):
    s = get_scalar(B1, pol)
    return eval_one(rows, s, stage2(B1, m, pol), t1_cost(B1, pol), mode)


def _eval_N(args):
    N, mode = args
    rows = load(N)
    res = []
    for B1, pol, m in candidates(B1_GRID):
        t, fr, nc = ev(rows, B1, pol, m, mode)
        res.append((t, fr, nc, B1, pol, m))
    res.sort(key=lambda r: r[0])
    return N, res


LEGACY_DISPATCH = {N: (27 if N <= 20 else 47 if N <= 22 else 70 if N <= 24 else 85 if N <= 26
                       else 125 if N <= 29 else 165 if N <= 31 else 205) for N in range(10, 40)}


def evaluate(N0, N1, mode, out=None):
    lines = []
    with Pool() as pool:
        allres = pool.map(_eval_N, [(N, mode) for N in range(N0, N1 + 1)])
    for N, res in allres:
        best = res[0]
        lines.append('N %d (%s): best %.0f ns  B1 %d pol %s B2 %dxB1  found %.4f  curves %.2f'
                     % (N, mode, best[0], best[3], best[4], best[5], best[1], best[2]))
        for r in res[1:4]:
            lines.append('      %.0f ns (%.3f)  B1 %d pol %s B2 %dxB1  found %.4f  curves %.2f'
                         % (r[0], r[0] / best[0], r[3], r[4], r[5], r[1], r[2]))
        lb = LEGACY_DISPATCH[N]
        lg = [r for r in res if r[3] == lb and r[4] == 'lvl' and r[5] == legacy_mult(lb)][0]
        lines.append('   legacy dispatch B1 %d (x%d): %.0f ns (%.3f)  found %.4f  curves %.2f'
                     % (lb, lg[5], lg[0], lg[0] / best[0], lg[1], lg[2]))
        lv = sorted([r for r in res if r[4] == 'lvl' and r[5] == legacy_mult(r[3])],
                    key=lambda r: r[0])
        lines.append('   best legacy level B1 %d: %.3f' % (lv[0][3], lv[0][0] / best[0]))
        byb = {}
        for r in res:
            if r[3] not in byb:
                byb[r[3]] = r
        lines.append('   best per B1: ' + '  '.join('%d:%.2f(x%d)' % (b, byb[b][0] / best[0], byb[b][5])
                                                for b in sorted(byb)))
    for l in lines:
        print(l)
    if out:
        with open(out, 'w') as f:
            f.write('\n'.join(lines) + '\n')


def _decomp_N(args):
    N, mode = args
    rows = load(N)
    cache = {}

    def e(B1, pol, m):
        if (B1, pol, m) not in cache:
            cache[(B1, pol, m)] = ev(rows, B1, pol, m, mode)[0]
        return cache[(B1, pol, m)]

    lb = LEGACY_DISPATCH[N]
    legacy = e(lb, 'lvl', legacy_mult(lb))
    k0 = lambda x: x[0]     # noqa: E731
    lm = min([(e(lb, 'lvl', m), m) for m in MULTS], key=k0)
    lpm = min([(e(lb, pol, m), pol, m) for pol in POLS + ['lvl'] for m in MULTS], key=k0)
    anyl = min([(e(B1, pol, m), B1, pol, m) for B1, pol, m in candidates(LEVELS)], key=k0)
    anylv = min([(e(B1, 'lvl', legacy_mult(B1)), B1) for B1 in LEVELS], key=k0)
    best = min([(e(B1, pol, m), B1, pol, m) for B1, pol, m in candidates(B1_GRID)], key=k0)
    return ('N %d (%s): legacy dispatch B1 %d %.0f | best legacy level %d %.3f | +mult x%d %.3f '
            '| +pol %s x%d %.3f | legacy B1s free %d %s x%d %.3f | best %d %s x%d %.3f'
            % (N, mode, lb, legacy, anylv[1], legacy / anylv[0], lm[1], legacy / lm[0],
               lpm[1], lpm[2], legacy / lpm[0], anyl[1], anyl[2], anyl[3], legacy / anyl[0],
               best[1], best[2], best[3], legacy / best[0]))


def decomp(N0, N1, mode, out=None):
    with Pool() as pool:
        lines = pool.map(_decomp_N, [(N, mode) for N in range(N0, N1 + 1)])
    for l in lines:
        print(l)
    if out:
        with open(out, 'w') as f:
            f.write('\n'.join(lines) + '\n')


def _plan_tab(args):
    N, mode, levels = args
    rows = load(N)
    return N, {(B1, pol, m): ev(rows, B1, pol, m, mode)[0] for B1, pol, m in candidates(levels)}


def plan(N0, N1, mode, levels, out=None):
    """one (policy, multiplier) per level and a level per N: pick per level
    the (pol, mult) with the least summed relative time over the N's where
    the level is (then) best, iterate; compare with the legacy dispatch."""
    Ns = list(range(N0, N1 + 1))
    with Pool() as pool:
        res = pool.map(_plan_tab, [(N, mode, levels) for N in Ns])
    tab = {N: t for N, t in res}
    bestN = {N: min(v for k, v in tab[N].items()) for N in Ns}
    opts = [(pol, m) for pol in POLS + ['lvl'] for m in MULTS]
    cfg = {B1: ('lvl' if B1 in LEVELS else (4, 1), legacy_mult(B1)) for B1 in levels}

    def tt(N, B1):
        return tab[N].get((B1,) + cfg[B1], 1e30)

    for it in range(5):
        disp = {N: min(levels, key=lambda b: tt(N, b)) for N in Ns}
        for B1 in levels:
            served = [N for N in Ns if disp[N] == B1] or \
                [min(Ns, key=lambda N: tt(N, B1) / bestN[N])]
            cfg[B1] = min((o for o in opts if (B1,) + o in tab[Ns[0]]),
                          key=lambda o: sum(tab[N][(B1,) + o] / bestN[N] for N in served))
    disp = {N: min(levels, key=lambda b: tt(N, b)) for N in Ns}
    lines = ['mode %s levels: ' % mode + '  '.join('%d: pol %s x%d' % (b, cfg[b][0], cfg[b][1])
                                                for b in levels)]
    for N in Ns:
        lb = LEGACY_DISPATCH[N]
        lt = tab[N].get((lb, 'lvl', legacy_mult(lb)))
        nt = tt(N, disp[N])
        lines.append('N %d: legacy B1 %d %.0f ns | new B1 %d %.0f ns  gain %.3f | best-possible %.3f'
                     % (N, lb, lt if lt else 0, disp[N], nt, (lt or nt) / nt, (lt or nt) / bestN[N]))
    for l in lines:
        print(l)
    if out:
        with open(out, 'w') as f:
            f.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'gen':
        gen(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]),
            int(sys.argv[6]) if len(sys.argv) > 6 else 1)
    elif cmd == 'eval':
        evaluate(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4],
                 sys.argv[5] if len(sys.argv) > 5 else None)
    elif cmd == 'decomp':
        decomp(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4],
               sys.argv[5] if len(sys.argv) > 5 else None)
    elif cmd == 'plan':
        plan(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4],
             [int(v) for v in sys.argv[5].split(',')], sys.argv[6] if len(sys.argv) > 6 else None)
    else:
        print(__doc__)
