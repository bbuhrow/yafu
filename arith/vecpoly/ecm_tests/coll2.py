# coll2.py p sigma S D dF blocks M: list the zero denominators mod ord(Q) in
# yafu's Brent-Suyama tables (1-D and 2-D setups, main steps): each one is a
# (genuine) stage-2 find of p that the oracle's f(jD) = +-f(k) test misses.
import sys, math
sys.path.insert(0, '.')
import bs_oracle as orc
p, sigma, S, D, dF, blocks, M = [int(x) for x in sys.argv[1:8]]
B1 = int(sys.argv[8]) if len(sys.argv) > 8 else 3000
x0, A, g = orc.suyama(sigma, p)
r0 = orc.point_order((x0, 1), A, g, p)
r = r0 // math.gcd(r0, orc.stage1_exponent(B1))
print('ord(Q) =', r)
f = orc.make_f(S); deg = abs(S)
def fd(i0, delta, d):
    v = [f(i0 + t * delta) for t in range(d + 1)]
    for e in range(1, d + 1):
        for t in range(d, e - 1, -1):
            v[t] -= v[t - 1]
    return v
def bad(a, b):
    return a != b and ((a - b) % r == 0 or (a + b) % r == 0)
hits = []
def pts(tag, vals):
    for v in vals:
        if v % r == 0: hits.append((tag, 'ladder T=0'))
def steps(tag, rows, n):
    for m in range(n):
        for row in rows:
            for e in range(len(row) - 1):
                if bad(row[e], row[e + 1]): hits.append((tag, 'step', m))
                if row[e] == row[e + 1] and (2 * row[e]) % r == 0: hits.append((tag, 'dbl', m))
            for e in range(len(row) - 1):
                row[e] += row[e + 1]
def family(tag, base, rstep, C, delta, twod, nsteps):
    cols = []
    if twod:
        g_ = [fd(base + k * rstep, delta, deg - k) for k in range(deg + 1)]
        tri = []
        for e in range(deg + 1):
            row = [g_[k][e] for k in range(deg - e + 1)]
            for j in range(1, len(row)):
                for k in range(len(row) - 1, j - 1, -1):
                    row[k] -= row[k - 1]
            tri.append(row)
        pts(tag + '-tri', [v for row in tri for v in row])
        for c in range(C):
            cols.append([row[0] for row in tri])
            if c + 1 < C: steps(tag + '-tri', tri, 1)
    else:
        for c in range(C):
            cols.append(fd(base + c * rstep, delta, deg))
        pts(tag, [v for col in cols for v in col])
    steps(tag, cols, nsteps)
for twod in (0, 1):
    hits = []
    family('giant', ((B1 + D // 2) // D) * D, D, M, M * D, twod, blocks * dF // M)
    for b in (1, 5):
        family('baby%d' % b, b, 6, M // 2, 3 * M, twod, (D // 2) // (3 * M) + 1)
    print('2-D' if twod else '1-D', hits[:10], len(hits))
