"""
Exact oracle for AVX-ECM stage 2 (plain and Brent-Suyama), small p.

For each curve (Suyama parametrization, sigma as yafu assigns them) it
computes the order of the starting point modulo p by baby-step giant-step,
applies yafu's stage-1 exponent, and decides whether the FFT stage 2 with
root polynomial f must find p: with Q the stage-1 point of order r, stage 2
finds p iff f(j D) = +-f(k) (mod r) for some giant j and baby k (the
product of x([f(jD)]Q) - x([f(k)]Q) vanishes mod p).

  python3 bs_oracle.py selftest
  python3 bs_oracle.py lanes p B1 D dF blocks S sigma0 nbatch
      S: 1 plain, S > 1 power x^S, S < -1 Dickson D_|S|(x, -1)
      prints, per batch of 8 curves (sigma0 + 8 b + lane), the lanes found
      in stage 1 and in stage 2
"""
import math
import random
import sys


# ---- Montgomery curve g y^2 = x^3 + A x^2 + x over F_p, affine ----
def m_add(P, Q, A, g, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + 2 * A * x1 + 1) * pow(2 * g * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (g * lam * lam - A - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def m_neg(P, p):
    return None if P is None else (P[0], (-P[1]) % p)


def m_mul(k, P, A, g, p):
    if k < 0:
        return m_mul(-k, m_neg(P, p), A, g, p)
    R = None
    while k:
        if k & 1:
            R = m_add(R, P, A, g, p)
        P = m_add(P, P, A, g, p)
        k >>= 1
    return R


# x-only ladder on (X:Z), as the C code: returns ([k]P, [k+1]P) projective
def x_ladder(k, x0, A, p):
    a24 = (A + 2) * pow(4, -1, p) % p

    def dbl(X, Z):
        s, d = (X + Z) ** 2 % p, (X - Z) ** 2 % p
        t = (s - d) % p
        return s * d % p, t * (d + a24 * t) % p

    def add(X1, Z1, X2, Z2, Xd, Zd):
        u = (X1 - Z1) * (X2 + Z2) % p
        v = (X1 + Z1) * (X2 - Z2) % p
        return Zd * (u + v) ** 2 % p, Xd * (u - v) ** 2 % p

    R0 = (x0, 1)
    R1 = dbl(x0, 1)
    for b in bin(k)[3:]:
        if b == '1':
            R0 = add(*R0, *R1, x0, 1)
            R1 = dbl(*R1)
        else:
            R1 = add(*R0, *R1, x0, 1)
            R0 = dbl(*R0)
    return R0, R1


# Okeya-Sakurai y-recovery, projective form used in the C code:
# P = (x0, 1) on g y^2 = ..., (X1:Z1) = [k]P, (X2:Z2) = [k+1]P.
# Returns (xn, yn, W) with [k]P = (xn / W, yn / W).
def y_recover(x0, A, g, X1, Z1, X2, Z2, p):
    N = (Z2 * ((x0 * X1 + Z1) * (x0 * Z1 + X1 + 2 * A * Z1) - 2 * A * Z1 * Z1)
         - (x0 * Z1 - X1) ** 2 * X2) % p
    W = 2 * g * Z1 * Z1 * Z2 % p
    xn = X1 * (2 * g * Z1 * Z2) % p
    return xn, N, W


# ---- Suyama curve from sigma, mod p ----
def suyama(sigma, p):
    u = (sigma * sigma - 5) % p
    v = 4 * sigma % p
    x0 = u ** 3 * pow(v ** 3, -1, p) % p
    A = ((v - u) ** 3 * (3 * u + v) * pow(4 * u ** 3 * v, -1, p) - 2) % p
    g = (x0 ** 3 + A * x0 * x0 + x0) % p
    return x0, A, g


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def rho(n):
    if n % 2 == 0:
        return 2
    while True:
        c = random.randrange(1, n)
        f = lambda x: (x * x + c) % n
        x = y = random.randrange(2, n)
        d = 1
        while d == 1:
            x = f(x)
            y = f(f(y))
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d


def factor(n, out=None):
    out = {} if out is None else out
    if n == 1:
        return out
    if is_prime(n):
        out[n] = out.get(n, 0) + 1
        return out
    for q in range(2, 1000):
        if n % q == 0:
            out[q] = out.get(q, 0) + 1
            return factor(n // q, out)
    d = rho(n)
    factor(d, out)
    factor(n // d, out)
    return out


def point_order(P, A, g, p):
    """order of P on the curve (Hasse interval, BSGS)"""
    w = 4 * math.isqrt(p) + 4
    L = p + 1 - w // 2
    m = math.isqrt(w) + 1
    baby = {}
    R = None
    for i in range(m):
        if R is not None:
            baby.setdefault(R, i)
        else:
            baby.setdefault('O', i)
        R = m_add(R, P, A, g, p)
    step = m_mul(m, P, A, g, p)
    G = m_mul(L, P, A, g, p)
    M = None
    for t in range(m + 2):
        key = 'O' if G is None else m_neg(G, p)
        if key in baby:
            M = L + t * m + baby[key]
            break
        G = m_add(G, step, A, g, p)
    assert M is not None and m_mul(M, P, A, g, p) is None
    for q in factor(M):
        while M % q == 0 and m_mul(M // q, P, A, g, p) is None:
            M //= q
    return M


def stage1_exponent(B1):
    k = 1
    q = 2
    while q < B1:
        k *= 2
        q *= 2
    for p in range(3, B1):
        if is_prime(p):
            c = p
            k *= p
            while c * p < B1:
                c *= p
                k *= p
    return k


def make_f(S):
    if abs(S) <= 1:
        return lambda x: x
    if S > 1:
        return lambda x: x ** S
    s = -S

    def dickson(x):
        d0, d1 = 2, x
        for _ in range(s - 1):
            d0, d1 = d1, x * d1 + d0          # D_n = x D_{n-1} - a D_{n-2}, a = -1
        return d1
    return dickson


def covered(r, f, D, j0, ngiant):
    babies = set()
    for k in range(1, D // 2):
        if math.gcd(k, D) == 1:
            v = f(k) % r
            babies.add(v)
            babies.add((-v) % r)
    for j in range(j0, j0 + ngiant):
        if f(j * D) % r in babies:
            return True
    return False


def lanes(p, B1, D, dF, blocks, S, sigma0, nbatch):
    k1 = stage1_exponent(B1)
    f = make_f(S)
    j0 = max(1, (B1 + D // 2) // D)
    for b in range(nbatch):
        s1, s2 = [], []
        for lane in range(8):
            sigma = sigma0 + 8 * b + lane
            x0, A, g = suyama(sigma, p)
            r0 = point_order((x0, 1), A, g, p)
            r = r0 // math.gcd(r0, k1)
            if r == 1:
                s1.append(lane)
                s2.append(lane)     # Q is the identity mod p: stage 2's gcd has p too
            elif covered(r, f, D, j0, blocks * dF):
                s2.append(lane)
        print('batch %d sigma0 %d: stage1 %s stage2 %s' %
              (b, sigma0 + 8 * b, ' '.join(map(str, s1)), ' '.join(map(str, s2))))


def selftest():
    random.seed(1)
    p = 1000003
    while not is_prime(p):
        p += 2
    p = (1 << 40) + 15
    while not is_prime(p):
        p += 2
    for _ in range(50):
        sigma = random.randrange(6, 1 << 30)
        x0, A, g = suyama(sigma, p)
        P = (x0, 1)
        assert (g * 1 - (x0 ** 3 + A * x0 * x0 + x0)) % p == 0
        k = random.randrange(2, 1 << 60)
        (X1, Z1), (X2, Z2) = x_ladder(k, x0, A, p)
        R = m_mul(k, P, A, g, p)
        if R is None or Z1 % p == 0:
            continue
        # x-only ladder agrees with the affine group law on the g-curve
        assert X1 * pow(Z1, -1, p) % p == R[0]
        xn, yn, W = y_recover(x0, A, g, X1, Z1, X2, Z2, p)
        Wi = pow(W, -1, p)
        assert (xn * Wi % p, yn * Wi % p) == R, 'y-recovery'
    print('selftest: ladder vs affine law and y-recovery OK')


if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        selftest()
    elif sys.argv[1] == 'lanes':
        a = [int(x) for x in sys.argv[2:]]
        lanes(*a)
