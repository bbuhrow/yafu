# Exact reference: which Suyama sigmas find p with B1, B2 (prime continuation)
import sys, math

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def ladder(k, X, Z, a24, p):
    # x-only Montgomery ladder: returns [k](X:Z)
    if k == 1: return X, Z
    x0, z0 = X, Z
    # R0 = P, R1 = 2P
    def dbl(x, z):
        s = (x + z) % p; d = (x - z) % p
        ss = s * s % p; dd = d * d % p; t = (ss - dd) % p
        return ss * dd % p, t * (dd + a24 * t) % p
    def add(x1, z1, x2, z2, xd, zd):
        u = (x1 - z1) * (x2 + z2) % p; v = (x1 + z1) * (x2 - z2) % p
        return zd * (u + v) ** 2 % p, xd * (u - v) ** 2 % p
    r0 = (x0, z0); r1 = dbl(x0, z0)
    for bit in bin(k)[3:]:
        if bit == '1':
            r0 = add(*r0, *r1, x0, z0); r1 = dbl(*r1)
        else:
            r1 = add(*r0, *r1, x0, z0); r0 = dbl(*r0)
    return r0

def finds(sig, p, B1, B2, P1, P2):
    u = (sig * sig - 5) % p; v = 4 * sig % p
    X = pow(u, 3, p); Z = pow(v, 3, p)
    num = pow(v - u, 3, p) * (3 * u + v) % p
    den = 16 * pow(u, 3, p) * v % p
    if den == 0: return 'bad'
    a24 = num * pow(den, -1, p) % p
    for q in P1:
        qq = q
        while qq * q <= B1: qq *= q
        X, Z = ladder(qq, X, Z, a24, p)
    if Z % p == 0: return 's1'
    # baby steps [j]Q, giant steps [mD]Q; q = mD +- j found iff
    # X_m Z_j - X_j Z_m == 0 (mod p)
    D = 2310
    baby = {}
    for j in range(1, D // 2 + 1, 2):
        if math.gcd(j, D) == 1:
            baby[j] = ladder(j, X, Z, a24, p)
    G = ladder(D, X, Z, a24, p)
    def add(x1, z1, x2, z2, xd, zd):
        u = (x1 - z1) * (x2 + z2) % p; v = (x1 + z1) * (x2 - z2) % p
        return zd * (u + v) ** 2 % p, xd * (u - v) ** 2 % p
    mlo = (B1 + D // 2) // D; mhi = (B2 + D // 2) // D + 1
    giant = {}
    prev = ladder(max(mlo - 1, 1) * D, X, Z, a24, p) if mlo > 1 else None
    cur = ladder(max(mlo, 1) * D, X, Z, a24, p)
    m = max(mlo, 1)
    giant[m] = cur
    while m < mhi:
        nxt = add(*cur, *G, *prev) if prev is not None else ladder((m + 1) * D, X, Z, a24, p)
        prev, cur = cur, nxt; m += 1; giant[m] = cur
    acc = 1
    for q in P2:
        mm = (q + D // 2) // D; j = abs(q - mm * D)
        xm, zm = giant[mm]; xj, zj = baby[j]
        acc = acc * (xm * zj - xj * zm) % p
    return 's2' if acc == 0 else None

p = int(sys.argv[1]); B1 = int(sys.argv[2]); B2 = int(sys.argv[3])
s0 = int(sys.argv[4]); kmax = int(sys.argv[5])
ps = primes_upto(B2)
P1 = [q for q in ps if q <= B1]; P2 = [q for q in ps if q > B1]
row = []
for k in range(kmax + 1):
    hits = []
    for j in range(8):
        r = finds(s0 + 8 * k + j, p, B1, B2, P1, P2)
        if r: hits.append(f"{j}{r}")
    row.append(",".join(hits) if hits else ".")
print("ref:", " ".join(row))
