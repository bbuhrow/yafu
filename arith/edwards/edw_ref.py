"""
Phase 0 reference for AVX-ECM stage 1 on a = -1 twisted Edwards curves
(arith/avxecm_edwards_plan.md).  Pure Python, exact arithmetic mod p.

  python edw_ref.py selftest
      group law (extended a = -1 formulas vs the affine law), signed-window
      scalar multiplication, the Edwards -> Montgomery hand-off (x-only ladder
      with A24 = 1/(1+d)), curve builders, torsion
  python edw_ref.py lanes fam p B1 D dF blocks S m0 nbatch
      exact oracle, like ../vecpoly/ecm_tests/bs_oracle.py lanes: per batch of
      8 curves (m0 + 8 b + lane) the lanes found at curve build (stage 0), in
      stage 1 and in stage 2.  fam: suyama | iv | j1 | iv-mont (Phase 1:
      the Montgomery image through PRAC)
  python edw_ref.py yield bits B1 B2 ncurves [seed]
      paired yield comparison of the curve families on random primes

Curve families ("m" plays the role of sigma):
  suyama  yafu's Montgomery curves from sigma = m
  iv      Gelin-Kleinjung-Lenstra family (iv): [m]Q on E': y^2 = x^3 - x^2
          - 9x + 9, Q = (5, 8); t = (4x + 4)/(y - 4); Thm 3.5 j = 4
  j1      GKL Thm 3.5 j = 1 with t = m (generic Z/2 x Z/4, a = -1)
All Edwards curves are -x^2 + y^2 = 1 + d x^2 y^2 with d = -e^4.
"""
import math
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'vecpoly', 'ecm_tests'))
import bs_oracle as bo                      # noqa: E402


class Find(Exception):
    """a non-invertible denominator mod p: the curve build / chain finds p"""


def inv(v, p):
    v %= p
    if v == 0:
        raise Find()
    return pow(v, -1, p)


# ---------------------------------------------------------------------------
# twisted Edwards a x^2 + y^2 = 1 + d x^2 y^2, affine (None-free: identity (0,1))
# ---------------------------------------------------------------------------
def ed_add(P, Q, a, d, p):
    x1, y1 = P
    x2, y2 = Q
    k = d * x1 * x2 * y1 * y2 % p
    x3 = (x1 * y2 + y1 * x2) * inv(1 + k, p) % p
    y3 = (y1 * y2 - a * x1 * x2) * inv(1 - k, p) % p
    return (x3, y3)


def ed_mul(k, P, a, d, p):
    if k < 0:
        return ed_mul(-k, ((-P[0]) % p, P[1]), a, d, p)
    R = (0, 1)
    while k:
        if k & 1:
            R = ed_add(R, P, a, d, p)
        P = ed_add(P, P, a, d, p)
        k >>= 1
    return R


def on_curve(P, a, d, p):
    x, y = P
    return (a * x * x + y * y - 1 - d * x * x * y * y) % p == 0


# ---------------------------------------------------------------------------
# extended coordinates, a = -1 (Hisil et al.; EFD dbl-2008-hwcd,
# add-2008-hwcd-3).  Op counters: M, S (field), A (add/sub).
# ---------------------------------------------------------------------------
class Ops:
    def __init__(self):
        self.M = self.S = self.A = 0
        self.dbl = self.dblt = self.add = 0

    def __repr__(self):
        return 'M=%d S=%d A=%d (dbl %d, dbl_t %d, add %d)' % (
            self.M, self.S, self.A, self.dbl, self.dblt, self.add)


def x_dbl(P, p, ops, want_t):
    """[2]P, P = (X, Y, Z[, T]); 3M + 4S (+1M for T)"""
    X, Y, Z = P[0], P[1], P[2]
    A = X * X % p
    B = Y * Y % p
    C = 2 * Z * Z % p                       # S + A
    E = ((X + Y) ** 2 - A - B) % p          # S + 3A
    G = (B - A) % p                         # D + B with D = -A
    F = (G - C) % p
    H = (-A - B) % p                        # D - B
    ops.S += 4
    ops.A += 7
    ops.M += 3
    X3, Y3, Z3 = E * F % p, G * H % p, F * G % p
    if want_t:
        ops.M += 1
        ops.dblt += 1
        return (X3, Y3, Z3, E * H % p)
    ops.dbl += 1
    return (X3, Y3, Z3)


def x_cache(P, d, p, ops):
    """cached form (Y-X, Y+X, 2d T, 2 Z) of an extended point; 1M + 3A"""
    X, Y, Z, T = P
    ops.M += 1
    ops.A += 3
    return ((Y - X) % p, (Y + X) % p, 2 * d * T % p, 2 * Z % p)


def x_neg_cached(C, p):
    """-(X, Y, Z, T) = (-X, Y, Z, -T): swap Y-X / Y+X, negate 2dT (free in a
    vector implementation: the subtraction order changes)"""
    return (C[1], C[0], (-C[2]) % p, C[3])


def x_add(P, C, p, ops, want_t):
    """P + cached C; 8M with T3, 7M without"""
    X1, Y1, Z1, T1 = P
    A = (Y1 - X1) * C[0] % p
    B = (Y1 + X1) * C[1] % p
    Cc = T1 * C[2] % p
    D = Z1 * C[3] % p
    E, F, G, H = (B - A) % p, (D - Cc) % p, (D + Cc) % p, (B + A) % p
    ops.M += 7
    ops.A += 6
    ops.add += 1
    if want_t:
        ops.M += 1
        return (E * F % p, G * H % p, F * G % p, E * H % p)
    return (E * F % p, G * H % p, F * G % p)


def to_affine(P, p):
    zi = inv(P[2], p)
    return (P[0] * zi % p, P[1] * zi % p)


# ---------------------------------------------------------------------------
# signed windows (width-w NAF), shared by all lanes in the C version
# ---------------------------------------------------------------------------
def wnaf(k, w):
    """digits, least significant first, in {0, +-1, +-3, ..., +-(2^(w-1)-1)}"""
    out = []
    half, full = 1 << (w - 1), 1 << w
    while k:
        if k & 1:
            r = k & (full - 1)
            if r >= half:
                r -= full
            k -= r
            out.append(r)
        else:
            out.append(0)
        k >>= 1
    return out


def x_scalar(k, P, d, p, w, ops):
    """[k]P for P extended (X, Y, Z, T), k > 0, with a width-w NAF and a table
    of cached odd multiples P, 3P, ..., (2^(w-1)-1)P.  A doubling computes T
    only when an addition follows (Hisil's mixing); an addition never computes
    T (a doubling always follows).  Returns (X, Y, Z) or (X, Y, Z, T)."""
    digs = wnaf(k, w)
    ext = [P]                               # extended odd multiples (top digit)
    tab = [x_cache(P, d, p, ops)]
    if w > 2:
        P2 = x_dbl(P, p, ops, True)
        c2 = x_cache(P2, d, p, ops)
        cur = P
        for _ in range(1, 1 << (w - 2)):
            cur = x_add(cur, c2, p, ops, True)
            ext.append(cur)
            tab.append(x_cache(cur, d, p, ops))
    R = ext[digs[-1] >> 1]                  # the top digit is positive
    for i in range(len(digs) - 2, -1, -1):
        v = digs[i]
        R = x_dbl(R, p, ops, v != 0)
        if v:
            c = tab[abs(v) >> 1]
            R = x_add(R, c if v > 0 else x_neg_cached(c, p), p, ops, False)
    return R


def stage1_ops(k, w):
    """op counts of x_scalar for a scalar k (no arithmetic): dbl / dbl_t / add"""
    digs = wnaf(k, w)
    ntab = (1 << (w - 2)) if w > 2 else 1
    nadd = sum(1 for v in digs[:-1] if v)
    ndbl_t = nadd
    ndbl = len(digs) - 1 - ndbl_t
    return ndbl, ndbl_t, nadd, ntab


# ---------------------------------------------------------------------------
# curve builders, mod p.  Return (d, x, y) for -x^2 + y^2 = 1 + d x^2 y^2,
# raising Find when a denominator vanishes mod p (yafu: a factor at build).
# ---------------------------------------------------------------------------
def w_add(P, Q, a2, a4, p):
    """y^2 = x^3 + a2 x^2 + a4 x + a6, affine; None = infinity"""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) * inv(2 * y1, p) % p
    else:
        lam = (y2 - y1) * inv(x2 - x1, p) % p
    x3 = (lam * lam - a2 - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def w_mul(k, P, a2, a4, p):
    """[k]P, left to right like vec_build_one_curve_edwards (avxecm.c):
    every step inverts a denominator, so an intermediate point at infinity
    mod p, or an addition R + P with x(R) = x(P) mod p, is a Find"""
    R = P
    for b in bin(k)[3:]:
        if R[1] % p == 0:
            raise Find()                    # 2y = 0: doubling to infinity
        R = w_add(R, R, a2, a4, p)
        if b == '1':
            if (R[0] - P[0]) % p == 0:
                raise Find()
            R = w_add(R, P, a2, a4, p)
    return R


def curve_j4_t(t, p):
    t %= p
    e = (t * t + 4 * t) * inv(t * t - 4, p) % p
    x1 = (2 * t ** 3 + 2 * t * t - 8 * t - 8) * inv(t ** 4 + 6 * t ** 3 + 12 * t * t + 16 * t, p) % p
    y1 = ((t ** 6 + 6 * t ** 5 + 10 * t ** 4 - 16 * t ** 3 - 48 * t * t - 32 * t - 32)
          * inv(t ** 6 + 6 * t ** 5 + 10 * t ** 4 + 16 * t ** 3 + 48 * t * t + 64 * t, p) % p)
    return e, x1, y1


def curve_iv(m, p):
    Q = w_mul(m, (5, 8), -1, -9, p)
    if Q is None:
        raise Find()
    x, y = Q
    t = (4 * x + 4) * inv(y - 4, p) % p
    e, x1, y1 = curve_j4_t(t, p)
    d = (-pow(e, 4, p)) % p
    return d, x1, y1


def curve_j1(m, p):
    t = m % p
    e = 3 * (t * t - 1) * inv(8 * t, p) % p
    x1 = 128 * t ** 3 * inv(27 * t ** 6 + 63 * t ** 4 - 63 * t * t - 27, p) % p
    y1 = (9 * t ** 4 - 2 * t * t + 9) * inv(9 * t ** 4 - 9, p) % p
    d = (-pow(e, 4, p)) % p
    return d, x1, y1


CURVES = {'iv': curve_iv, 'j1': curve_j1}


def edw_to_mont(d, x, y, p, a=-1):
    """Montgomery image: (A, B, u, v) with B v^2 = u^3 + A u^2 + u"""
    A = 2 * (a + d) * inv(a - d, p) % p
    B = 4 * inv(a - d, p) % p
    u = (1 + y) * inv(1 - y, p) % p
    v = u * inv(x, p) % p
    return A, B, u, v


def curve_point(fam, m, p):
    """(A, g, (x, y)) of the starting point's Montgomery image, for the order
    computation with bs_oracle (g y^2 = x^3 + A x^2 + x).  Raises Find."""
    if fam == 'suyama':
        x0, A, g = bo.suyama(m, p)
        return A, g, (x0, 1)
    d, x, y = CURVES[fam](m, p)
    if (1 + d) % p == 0:
        raise Find()                                # A24 = 1/(1+d) needed
    A, B, u, v = edw_to_mont(d, x, y, p)
    return A, B, (u, v)


# ---------------------------------------------------------------------------
# oracle
# ---------------------------------------------------------------------------
def stage1_r(fam, m, p, k1):
    """order of the stage-1 point, or 0 when the build finds p.  For Edwards a
    stage-1 point of order 2 is (0, -1): X = 0, found by the stage-1 gcd."""
    try:
        A, g, P = curve_point(fam, m, p)
    except Find:
        return 0, None
    r0 = bo.point_order(P, A, g, p)
    return r0 // math.gcd(r0, k1), (A, g, P)


def prac_stage1_z(x0, A, p, B1):
    """Z mod p after yafu's x-only stage 1 on (x0 : 1): doublings for 2^k <
    B1, then vec_prac(q) for every prime power q^j < B1 (ORIG_PRAC off), with
    the same chains (edw_model replays vec_prac's choices).  It differs from a
    ladder by s when the point's remaining order r is small: a PRAC chain for
    a later prime q > r can meet a difference that is the identity mod p, and
    the differential addition then collapses to Z = 0 -- a stage-1 find that
    the order alone does not predict."""
    import edw_model as em
    a24 = (A + 2) * pow(4, -1, p) % p

    def dup(P):
        s, d = (P[0] + P[1]) ** 2 % p, (P[0] - P[1]) ** 2 % p
        t = (s - d) % p
        return s * d % p, t * (d + a24 * t) % p

    def add(P1, P2, Pd):
        u = (P1[0] - P1[1]) * (P2[0] + P2[1]) % p
        v = (P1[0] + P1[1]) * (P2[0] - P2[1]) % p
        return Pd[1] * (u + v) ** 2 % p, Pd[0] * (u - v) ** 2 % p

    def prac(P, c):
        best, cmin = 0, 1e300
        for i, v in enumerate(em.VAL):
            cs = em.lucas_cost(c, v)
            if cs < cmin:
                cmin, best = cs, i
        r = int(c * em.VAL[best] + 0.5)
        d, e = c - r, 2 * r - c
        A_, B_, C_ = dup(P), P, P
        while d != e:
            if d < e:
                d, e = e, d
                A_, B_ = B_, A_
            if (d + 3) // 4 <= e:
                d -= e
                B_, C_ = add(B_, A_, C_), B_
            elif (d + e) % 2 == 0:
                d = (d - e) // 2
                B_ = add(B_, A_, C_)
                A_ = dup(A_)
            elif d % 2 == 0:
                d //= 2
                C_ = add(C_, A_, B_)
                A_ = dup(A_)
            else:
                e //= 2
                C_ = add(C_, B_, A_)
                B_ = dup(B_)
        return add(A_, B_, C_)

    P = (x0 % p, 1)
    q = 2
    while q < B1:
        P = dup(P)
        q *= 2
    for q in em.primes_below(B1)[1:]:
        c = 1
        while True:
            P = prac(P, q)
            c *= q
            if c * q >= B1:
                break
    return P[1] % p


def lanes(fam, p, B1, D, dF, blocks, S, m0, nbatch):
    """fam 'iv-mont': the curve runs as its Montgomery image through yafu's
    x-only PRAC stage 1 (Phase 1, -avxecm_curve 1).  For x-only stage 1
    (suyama, iv-mont) a remaining order 1 < r < B1 is checked by replaying
    the PRAC chains (prac_stage1_z)."""
    mont = fam.endswith('-mont')
    fam = fam.replace('-mont', '')
    xonly = mont or fam == 'suyama'
    k1 = bo.stage1_exponent(B1)
    f = bo.make_f(S)
    j0 = max(1, (B1 + D // 2) // D)
    for b in range(nbatch):
        s0, s1, s2 = [], [], []
        for lane in range(8):
            m = m0 + 8 * b + lane
            r, info = stage1_r(fam, m, p, k1)
            if r == 0:
                s0.append(lane)
                s1.append(lane)
                s2.append(lane)
            elif r == 1 or (r == 2 and not xonly):
                s1.append(lane)
                s2.append(lane)
            elif xonly and r < B1 and prac_stage1_z(info[2][0], info[0], p, B1) == 0:
                s1.append(lane)
                s2.append(lane)
            elif bo.covered(r, f, D, j0, blocks * dF):
                s2.append(lane)
        print('batch %d m0 %d: stage0 %s stage1 %s stage2 %s' % (
            b, m0 + 8 * b, ' '.join(map(str, s0)), ' '.join(map(str, s1)),
            ' '.join(map(str, s2))))


# ---------------------------------------------------------------------------
# yield comparison
# ---------------------------------------------------------------------------
def v_q(n, q):
    c = 0
    while n % q == 0:
        n //= q
        c += 1
    return c


def _yield_one(args):
    p, B1, B2, k1, ms = args
    out = {}
    for fam, m in ms.items():
        try:
            A, g, P = curve_point(fam, m, p)
        except Find:
            out[fam] = (1, 1, 0, 0)
            continue
        r0 = bo.point_order(P, A, g, p)
        r = r0 // math.gcd(r0, k1)
        s1 = int(r <= 2 if fam != 'suyama' else r == 1)
        s2 = int(s1 or (bo.is_prime(r) and r <= B2))
        out[fam] = (s1, s2, v_q(r0, 2), v_q(r0, 3))
    return out


def yield_cmp(bits, B1, B2, ncurves, seed=1):
    import multiprocessing as mp
    rng = random.Random(seed)
    k1 = bo.stage1_exponent(B1)
    jobs = []
    for _ in range(ncurves):
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        while not bo.is_prime(p):
            p += 2
        ms = {'suyama': rng.randrange(6, 1 << 31), 'iv': rng.randrange(1, 1 << 31),
              'j1': rng.randrange(5, 1 << 31)}
        jobs.append((p, B1, B2, k1, ms))
    with mp.Pool() as pool:
        res = pool.map(_yield_one, jobs, chunksize=64)
    print('bits %d B1 %d B2 %d curves %d (one random curve of each family per prime)'
          % (bits, B1, B2, ncurves))
    print('%-7s %9s %9s %9s %9s %9s' % ('family', 'stage1', 'stage1+2', 'ratio', 'E v2', 'E v3'))
    base = None
    for fam in ('suyama', 'iv', 'j1'):
        s1 = sum(r[fam][0] for r in res)
        s2 = sum(r[fam][1] for r in res)
        v2 = sum(r[fam][2] for r in res) / len(res)
        v3 = sum(r[fam][3] for r in res) / len(res)
        base = base or s2
        print('%-7s %9d %9d %9.4f %9.3f %9.3f' % (fam, s1, s2, s2 / base, v2, v3))
    # paired difference iv - suyama (std error of the count difference)
    dif = [r['iv'][1] - r['suyama'][1] for r in res]
    mu = sum(dif) / len(dif)
    sd = math.sqrt(sum((x - mu) ** 2 for x in dif) / (len(dif) - 1))
    print('iv - suyama finds: %+d (std error %.0f)' % (sum(dif), sd * math.sqrt(len(dif))))


# ---------------------------------------------------------------------------
# self test
# ---------------------------------------------------------------------------
def selftest():
    rng = random.Random(7)
    nchk = 0
    for bits in (40, 61):
        for _ in range(40):
            p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
            while not bo.is_prime(p):
                p += 2
            fam = rng.choice(['iv', 'j1'])
            m = rng.randrange(1, 1 << 20)
            try:
                d, x, y = CURVES[fam](m, p)
            except Find:
                continue
            P = (x, y)
            assert on_curve(P, -1, d, p), 'builder: point not on curve'
            # extended formulas vs the affine law
            Pe = (x, y, 1, x * y % p)
            for k in (2, 3, 5, 7, 12, 1000, rng.randrange(1, 1 << 64)):
                for w in (2, 3, 4, 6):
                    try:
                        R = x_scalar(k, Pe, d, p, w, Ops())
                        Ra = ed_mul(k, P, -1, d, p)
                        Rx = to_affine(R, p)
                    except Find:
                        continue
                    assert Rx == Ra, 'x_scalar k=%d w=%d' % (k, w)
                    assert len(R) == 3 or (R[3] * R[2] - R[0] * R[1]) % p == 0, "T invariant"
                    nchk += 1
            # hand-off: Montgomery x-only ladder with A24 = 1/(1+d) on
            # (Z+Y : Z-Y) matches the Edwards multiple
            A, B, u, v = edw_to_mont(d, x, y, p)
            a24 = inv(1 + d, p)
            assert (A + 2) * inv(4, p) % p == a24, 'A24 = 1/(1+d)'
            assert (B * v * v - (u ** 3 + A * u * u + u)) % p == 0, 'Montgomery image'
            k = rng.randrange(1, 1 << 40)
            R = x_scalar(k, Pe, d, p, 5, Ops())
            Xm, Zm = (R[2] + R[1]) % p, (R[2] - R[1]) % p
            (X1, Z1), _ = bo.x_ladder(k, u, A, p)
            assert (Xm * Z1 - X1 * Zm) % p == 0, 'Edwards -> Montgomery hand-off'
            # 2- and 4-torsion: (0,-1) order 2, (1/e, 1/e)-type points order 4
            assert ed_mul(2, (0, p - 1), -1, d, p) == (0, 1)
            nchk += 1
    # torsion: family (iv) has 16 | #E(F_p) (counted directly on small p),
    # j1 has 8 | #E
    for fam, need in (('iv', 16), ('j1', 8)):
        for _ in range(30):
            p = rng.randrange(1 << 11, 1 << 13) | 1
            while not bo.is_prime(p):
                p += 2
            try:
                d, x, y = CURVES[fam](rng.randrange(1, 1000), p)
            except Find:
                continue
            if d % p == 0 or (d + 1) % p == 0:
                continue
            n = count_points(d, p)
            assert n % need == 0, '%s: #E = %d not divisible by %d (p=%d)' % (fam, n, need, p)
            nchk += 1
    print('selftest OK (%d checks)' % nchk)


def count_points(d, p):
    """#E(F_p) of -x^2 + y^2 = 1 + d x^2 y^2 via its Montgomery model
    (projective, so points at infinity of the Edwards model are included)"""
    A = 2 * (d - 1) * inv(-1 - d, p) % p
    B = 4 * inv(-1 - d, p) % p
    n = 1                                   # infinity
    for u in range(p):
        r = (u ** 3 + A * u * u + u) * inv(B, p) % p
        n += 1 + legendre(r, p)
    return n


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'selftest':
        selftest()
    elif cmd == 'lanes':
        fam = sys.argv[2]
        lanes(fam, *[int(v) for v in sys.argv[3:]])
    elif cmd == 'yield':
        a = [int(v) for v in sys.argv[2:]]
        yield_cmp(*a)
