"""
Brent-Suyama extension: success probability and time per factor.

A stage 2 with root polynomial f (x^S or Dickson D_S(x, -1)) finds a curve
whose order is B1-smooth apart from one prime q when f(jD) = +-f(k) (mod q)
for a giant j and a baby k.  For q <= B2 that always happens (f(jD) - f(k)
contains jD - k and jD + k).  For q > B2 the expected number of pairs is
    E(q) = pairs * rbar_f / q,   pairs = nbaby * ngiant,
with rbar_f the average number of non-trivial x (mod q) with
f(x) = +-f(y), and the chance is 1 - exp(-E(q)).

  python3 bs_model.py rbar                     rbar_f for the polynomials
  python3 bs_model.py rates p B1 D dF blocks S ncurves
        oracle-exact stage-2 find rates (plain and S) on ncurves curves vs
        the model's prediction for this q distribution
  python3 bs_model.py table [--isa ...]        time per factor, plain vs BS
"""
import math
import random
import sys

import bs_oracle as orc
import s2model as sm


# ---- rbar_f by direct counting over random primes ----
def _primes_between(lo, hi, cnt, rng):
    out = []
    while len(out) < cnt:
        q = rng.randrange(lo, hi) | 1
        if orc.is_prime(q):
            out.append(q)
    return out


_RSAMP = {}


def rsamples(S, nq=150, ny=8, lo=3000, hi=9000, seed=7):
    """samples of the count of non-trivial x with f(x) = +-f(y) mod q, over
    random primes q and y (the count depends on q's residue class, so the
    chance 1 - exp(-pairs r / q) is averaged over the samples, not taken at
    the mean)"""
    if abs(S) <= 1:
        return [0]
    if S in _RSAMP:
        return _RSAMP[S]
    rng = random.Random(seed)
    f = orc.make_f(S)
    out = []
    for q in _primes_between(lo, hi, nq, rng):
        vals = [f(x) % q for x in range(q)]
        cnt = {}
        for v in vals:
            cnt[v] = cnt.get(v, 0) + 1
        for _ in range(ny):
            y = rng.randrange(1, q)
            fy = vals[y]
            c = cnt.get(fy, 0) + (cnt.get((-fy) % q, 0) if fy != 0 else 0)
            out.append(c - 2)     # x = y and x = -y are the trivial ones
    _RSAMP[S] = out
    return out


def rbar(S):
    r = rsamples(S)
    return sum(r) / len(r)


_RHIST = {}


def p_catch(S, E1):
    """mean of 1 - exp(-E1 r) over the r samples (E1 = pairs / q)"""
    if S not in _RHIST:
        r = rsamples(S)
        h = {}
        for x in r:
            h[x] = h.get(x, 0) + 1
        _RHIST[S] = [(x, c / len(r)) for x, c in h.items()]
    return sum(w * (1.0 - math.exp(-E1 * x)) for x, w in _RHIST[S])


# ---- probability with the extension ----
def prob_bs(pdigits, B1, B2, pairs, S, torsion=23.0, steps=200):
    """P(found) for a factor of pdigits digits, root polynomial S"""
    L = pdigits * math.log(10) - math.log(torsion)
    a, b = math.log(B1), math.log(B2)
    p = sm.rho(L / a)
    # one prime q in (B1, B2]
    h = (b - a) / steps
    s = 0.0
    for i in range(steps + 1):
        l = a + i * h
        w = 0.5 if i in (0, steps) else 1.0
        s += w * sm.rho((L - l) / a) / l
    p += s * h
    if abs(S) > 1 and L > b:
        # one prime q in (B2, x/B1 ...): caught with 1 - exp(-pairs rb / q)
        h = (L - b) / steps
        s = 0.0
        for i in range(steps + 1):
            l = b + i * h
            w = 0.5 if i in (0, steps) else 1.0
            s += w * sm.rho((L - l) / a) / l * p_catch(S, pairs / math.exp(l))
        p += s * h
    return p


# ---- oracle-exact rates on random curves ----
def rates(p, B1, D, dF, blocks, S, ncurves, seed=1):
    k1 = orc.stage1_exponent(B1)
    j0 = max(1, (B1 + D // 2) // D)
    ng = blocks * dF
    nb = sum(1 for k in range(1, D // 2) if math.gcd(k, D) == 1)
    f0, fs = orc.make_f(1), orc.make_f(S)
    B2 = (j0 + ng - 1) * D + D // 2 - 1
    rng = random.Random(seed)
    s1 = s2p = s2s = 0
    big = []                  # q > B2 from curves of prime order cofactor
    for _ in range(ncurves):
        sigma = rng.randrange(6, 1 << 40)
        x0, A, g = orc.suyama(sigma, p)
        try:
            r0 = orc.point_order((x0, 1), A, g, p)
        except (AssertionError, ValueError):
            continue
        r = r0 // math.gcd(r0, k1)
        if r == 1:
            s1 += 1
            continue
        cp = orc.covered(r, f0, D, j0, ng)
        cs = orc.covered(r, fs, D, j0, ng)
        s2p += cp
        s2s += cs
        if orc.is_prime(r) and r > B2:
            big.append((r, cs))
    rb = rbar(S)
    pred = sum(p_catch(S, nb * ng / q) for q, _ in big)
    got = sum(c for _, c in big)
    print('p %d, B1 %d, B2 %d (D %d, dF %d, %d blocks), f = %s, %d curves' %
          (p, B1, B2, D, dF, blocks, S, ncurves))
    print('  stage 1 finds %d; stage 2 finds plain %d, Brent-Suyama %d (+%.1f%%)' %
          (s1, s2p, s2s, 100.0 * (s2s - s2p) / max(1, s2p)))
    print('  curves with a single prime q > B2 left: %d; caught by BS %d, model predicts %.1f '
          '(rbar_f = %.2f, pairs = %d)' % (len(big), got, pred, rb, nb * ng))


# ---- time per factor: plain vs Brent-Suyama ----
def bs_m(t, blocks, deg):
    """s2_bs_m of avxecm_stage2fft.c"""
    dF, D = float(1 << t), sm.D_TAB[t]
    roots = dF * (blocks + 1.6)
    bits = deg * (math.log2((blocks * dF + 1.0) * D) + math.log2(D / 2.0))
    mopt = math.sqrt(roots * 8.0 * 12.0 / ((deg + 1.0) * 10.0 * bits))
    M = 2
    while M * 2 <= 64 and M * 2 <= dF and M * 1.41 < mopt:
        M *= 2
    return M


def cost_bs(t, n, k, lean, blocks, S, M=0):
    """s2_cost of avxecm_stage2fft.c with the Brent-Suyama terms"""
    dF = float(1 << t)
    T = sm.tcost(t, n) * dF
    rebuilt = sum(s % k for s in range(1, t)) / t
    groot, broot, init = 10.0 * dF, 17.0 * dF, 0.0
    deg = abs(S)
    if deg > 1:
        D = sm.D_TAB[t]
        M = M or bs_m(t, blocks, deg)
        # 7 mulmods + ~6 add/sub (O(n)) per affine addition; one batch
        # inversion (8 lanes x mpz_invert ~ 4 + 160/n mulmods) per M roots
        per = deg * (7.0 + 12.0 / n) + 1.0 + 8.0 * (4.0 + 160.0 / n) / M
        groot = dF * per
        broot = 1.6 * dF * per
        init = M * (deg + 1.0) * 10.0 * deg * (math.log2((blocks * dF + 1.0) * D) + math.log2(D / 2.0))
    return (blocks * (groot + T * (1.0 + (sm.MM_L if lean else sm.MM))) + init + broot
            + T * (1.0 + (sm.MI_L if lean else sm.MI) + (sm.EV_L if lean else sm.EV) + rebuilt))


def best(n, B1, pdig, cap, S):
    """least time per factor over configurations and block counts"""
    T1 = sm.t1_seconds(B1, n)
    out = None
    for t, k, lean, m in sm.configs(n, B1, cap):
        dF, D = 1 << t, sm.D_TAB[t]
        nb = sum(1 for kk in range(1, D // 2) if math.gcd(kk, D) == 1)
        bl = 1
        while True:
            T2 = sm.SCALE * cost_bs(t, n, k, lean, bl, S) * sm.mu_ns(n) * 1e-9
            if T2 > 4 * T1:
                break
            B2 = sm.covered_B2(t, B1, bl)
            e = (T1 + T2) / prob_bs(pdig, B1, B2, nb * bl * dF, S)
            if out is None or e < out[0]:
                out = (e, B2, T2, t, bl)
            bl = bl + 1 if bl < 8 else int(bl * 1.25)
    return T1, out


def table():
    polys = [1, 2, -3, 6, -6, -12, -30]
    print('rbar_f: ' + ', '.join('%s %.2f' % (name(S), rbar(S)) for S in polys[1:]))
    print('time per factor relative to the plain stage 2 (each with its best B2); ISA %s' % sm.ISA)
    print('bits  B1     p    MB  | ' + ' '.join('%8s' % name(S) for S in polys))
    for bits, B1, pd in [(400, 250000, 30), (400, 1000000, 35), (400, 3000000, 40),
                         (1024, 1000000, 35), (1024, 3000000, 40), (1024, 11000000, 45),
                         (2048, 3000000, 40), (2048, 11000000, 45)]:
        n = (bits + 3 + 51) // 52
        for cap in (64, 256):
            base = None
            row = []
            for S in polys:
                T1, (e, B2, T2, t, bl) = best(n, B1, pd, cap, S)
                if base is None:
                    base = e
                row.append('%8.3f' % (e / base))
            print('%4d %5.0e p%d %4d | %s' % (bits, B1, pd, cap, ' '.join(row)))


def name(S):
    return 'plain' if abs(S) <= 1 else ('x^%d' % S if S > 0 else 'D_%d' % -S)


if __name__ == '__main__':
    if sys.argv[1] == 'rbar':
        for S in (2, 3, 6, 12, -3, -6, -12, -30):
            print('%-6s rbar %.3f' % (name(S), rbar(S)))
    elif sys.argv[1] == 'rates':
        a = [int(x) for x in sys.argv[2:]]
        rates(*a)
    elif sys.argv[1] == 'table':
        table()
