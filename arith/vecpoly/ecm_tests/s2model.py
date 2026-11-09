"""
AVX-ECM FFT stage 2: memory / time model and the most efficient B2 for a
memory limit (arith/avxecm_stage2_fft_design.md section 11).

  python3 s2model.py                 tables for the design doc
  python3 s2model.py bits B1 MB      one case: every useful B2 at that limit

Memory follows avxecm_stage2fft.c (s2_mem); time is the s2_cost model
(template mulmods, fit to vp_test -b on IFMA) times the measured mulmod
time mu(n), scaled by 0.93 to match yafu runs.  Stage 1 time is
c1 * B1 * mu(n) (measured).  Success probability: Dickman rho for a group
order of size p / 23 (Suyama torsion), one prime in (B1, B2].
"""
import math
import sys

D_TAB = [0, 12, 30, 60, 120, 240, 510, 1020, 2310, 4620, 9240, 19110, 39270,
         79170, 158340, 324870]
GRP, SCHOOL_CAP, CHUNK = 4, 16, 256

# measured (Rocket Lake, IFMA, 1 thread): template mulmod ns and stage-1
# mulmods per unit of B1
MU_NS = {4: 17, 8: 59, 11: 110, 20: 352, 40: 1430}
C1 = {8: 16.9, 20: 13.9}


def mu_ns(n):
    if n in MU_NS:
        return MU_NS[n]
    return 352.0 * (n / 20.0) ** 2


def c1(n):
    return C1.get(n, 15.0)


def nprimes(n, t):
    return (104 * n + (t + 1) + 2) // 50 + 1


def levels(t, k):
    return (t + k - 1) // k


def mem_bytes(t, n, k, lean):
    """per thread, and the shared context"""
    K = nprimes(n, t + 1 - 1)
    dF = 1 << t
    L = 2 * dF
    cb = 64 * n
    work = (K + GRP) * L * 64 + (L // 2 + L // 4 + 4 * SCHOOL_CAP + 1) * cb + (6 * n + 16) * 64
    coef = (2 + levels(t, k) + (1 if lean else 0)) * dF * cb + (3 * min(dF, CHUNK) + 16) * cb
    modf = dF * cb + (0 if lean else K * 2 * dF * 64 + K * dF * 64)
    ctx = 2 * K * L * 8 + ((2 * K + 1) * n + (t + 2) * K * 2 + 3 * K) * 8
    return work + coef + modf, ctx


def tcost(t, n):
    return max(8.0, 8.0 * t - 42.8 + 4.05 * t * t / n)


def cost_mulmods(t, n, k, lean, blocks):
    dF = float(1 << t)
    T = tcost(t, n) * dF
    rebuilt = sum(s % k for s in range(1, t)) / t
    return (blocks * (10.0 * dF + T * (1.90 if lean else 1.78))
            + 17.0 * dF + T * (1.0 + (0.59 if lean else 0.82) + 1.88 + (0.1 if lean else 0) + rebuilt))


def t2_seconds(t, n, k, lean, blocks):
    return 0.93 * cost_mulmods(t, n, k, lean, blocks) * mu_ns(n) * 1e-9


def t1_seconds(B1, n):
    return c1(n) * B1 * mu_ns(n) * 1e-9


def covered_B2(t, B1, blocks):
    D, dF = D_TAB[t], 1 << t
    j0 = max(1, (B1 + D // 2) // D)
    return (j0 + blocks * dF - 1) * D + D // 2 - 1


# ---- Dickman rho ----
RHO_H = 1e-3
RHO_MAX = 30.0
_rho = None


def _rho_table():
    global _rho
    if _rho is not None:
        return _rho
    m = int(RHO_MAX / RHO_H) + 2
    r = [1.0] * m
    k1 = int(round(1.0 / RHO_H))
    for i in range(k1 + 1, m):
        u = i * RHO_H
        # rho'(u) = -rho(u-1)/u, trapezoid
        f0 = r[i - 1 - k1] / (u - RHO_H)
        f1 = r[i - k1] / u
        r[i] = r[i - 1] - 0.5 * RHO_H * (f0 + f1)
    _rho = r
    return r


def rho(u):
    if u <= 1.0:
        return 1.0
    r = _rho_table()
    x = u / RHO_H
    i = int(x)
    if i + 1 >= len(r):
        return 0.0
    f = x - i
    return r[i] * (1 - f) + r[i + 1] * f


def prob(pdigits, B1, B2, torsion=23.0):
    L = pdigits * math.log(10) - math.log(torsion)
    a, b = math.log(B1), math.log(B2)
    p = rho(L / a)
    steps = 200
    h = (b - a) / steps
    s = 0.0
    for i in range(steps + 1):
        l = a + i * h
        w = 0.5 if i in (0, steps) else 1.0
        s += w * rho((L - l) / a) / l
    return p + s * h


# ---- searches ----
def configs(n, B1, cap_mb, threads=1):
    cap = cap_mb * 1048576.0
    for t in range(3, 16):
        if B1 < D_TAB[t]:
            break
        for lean in (0, 1):
            for k in (1, 2, 3, 4):
                m, ctx = mem_bytes(t, n, k, lean)
                if m + ctx / threads <= cap:
                    yield t, k, lean, m


def best_for_B2(n, B1, B2, cap_mb):
    """what avxecm_s2fft_setup picks for a requested B2"""
    best = None
    for t, k, lean, m in configs(n, B1, cap_mb):
        dF, D = 1 << t, D_TAB[t]
        blocks = math.ceil((B2 - B1) / (dF * D))
        c = cost_mulmods(t, n, k, lean, blocks)
        if best is None or c < best[0]:
            best = (c, t, k, lean, blocks, m)
    return best


def efficient_points(n, B1, cap_mb, B2max):
    """(B2, T2, cfg) on the lower envelope of T2 versus covered B2"""
    pts = []
    for t, k, lean, m in configs(n, B1, cap_mb):
        dF, D = 1 << t, D_TAB[t]
        bl = 1
        while True:
            B2 = covered_B2(t, B1, bl)
            if B2 > B2max:
                break
            pts.append((B2, t2_seconds(t, n, k, lean, bl), (t, k, lean, bl, m)))
            bl = bl + 1 if bl < 8 else int(bl * 1.25)
    pts.sort()
    env = []
    for p in pts:
        while env and env[-1][1] >= p[1]:
            env.pop()
        env.append(p)
    return env


def optimum(n, B1, pdig, cap_mb):
    T1 = t1_seconds(B1, n)
    best = None
    for B2, T2, cfg in efficient_points(n, B1, cap_mb, B1 * 1e5):
        e = (T1 + T2) / prob(pdig, B1, B2)
        if best is None or e < best[0]:
            best = (e, B2, T2, cfg)
    return T1, best


def gmp_default_B2(B1):
    return (11.0 / 6.0 * B1) ** 1.43


def main():
    if len(sys.argv) == 4:
        bits, B1, cap = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
        n = (bits + 3 + 51) // 52
        T1 = t1_seconds(B1, n)
        print("n = %d, T1 = %.2f s; cheapest stage 2 for each covered B2 "
              "(B2 steps >= 15%%), with expected time per p35/p40 factor relative "
              "to the best" % (n, T1))
        env = [p for p in efficient_points(n, int(B1), cap, B1 * 1e5)
               if p[0] >= 10 * B1 and p[1] <= 4 * T1]
        eff = {pd: [(T1 + p[1]) / prob(pd, B1, p[0]) for p in env] for pd in (35, 40)}
        last = 0
        for i, (B2, T2, (t, k, lean, bl, m)) in enumerate(env):
            if B2 < 1.15 * last and i != len(env) - 1:
                continue
            last = B2
            print("B2 %9.3g  T2 %7.3f s  T2/T1 %5.2f  dF %5d  k=%d lean=%d blocks %4d  %6.1f MB"
                  "   p35 %5.3f  p40 %5.3f"
                  % (B2, T2, T2 / T1, 1 << t, k, lean, bl, m / 1048576,
                     eff[35][i] / min(eff[35]), eff[40][i] / min(eff[40])))
        return

    print("== memory per thread (MB), full / lean + every 3rd tree level ==")
    print("bits  " + "".join("   dF=%-6d" % (1 << t) for t in range(9, 16)))
    for bits in (200, 400, 700, 1024, 2048):
        n = (bits + 3 + 51) // 52
        row = []
        for t in range(9, 16):
            f, c = mem_bytes(t, n, 1, 0)
            l, _ = mem_bytes(t, n, 3, 1)
            row.append("%5.0f/%-5.0f" % ((f + c) / 1048576, (l + c) / 1048576))
        print("%4d  %s" % (bits, " ".join(row)))

    print()
    print("== most efficient B2 per memory limit (IFMA, 1 thread) ==")
    print("bits  B1     p   MB    | B2*      T2*/T1  dF    blocks | default B2  T2/T1  eff. loss | B2 within 2%")
    cases = [(400, 250000, 30), (400, 1000000, 35), (400, 3000000, 40),
             (1024, 250000, 30), (1024, 1000000, 35), (1024, 3000000, 40), (1024, 11000000, 45)]
    for bits, B1, pd in cases:
        n = (bits + 3 + 51) // 52
        for cap in (16, 64, 256, 1024):
            T1, (e, B2, T2, (t, k, lean, bl, m)) = optimum(n, B1, pd, cap)
            Bd = gmp_default_B2(B1)
            bd = best_for_B2(n, B1, Bd, cap)
            if bd is None:
                continue
            c, td, kd, ld, bld, md = bd
            Bdc = covered_B2(td, B1, bld)
            T2d = t2_seconds(td, n, kd, ld, bld)
            ed = (T1 + T2d) / prob(pd, B1, Bdc)
            near = [p[0] for p in efficient_points(n, B1, cap, B1 * 1e5)
                    if (T1 + p[1]) / prob(pd, B1, p[0]) <= 1.02 * e]
            print("%4d %4.0e p%d %5d | %8.2e  %5.2f  %5d  %4d   | %8.2e  %5.2f  %6.1f%%   | %.1e .. %.1e"
                  % (bits, B1, pd, cap, B2, T2 / T1, 1 << t, bl, Bdc, T2d / T1,
                     100.0 * (ed / e - 1.0), min(near), max(near)))


if __name__ == "__main__":
    main()
