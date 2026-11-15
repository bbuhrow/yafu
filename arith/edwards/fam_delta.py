"""
Extra smoothness of the ECM curve families (docs/ecm_prob.md).

GMP-ECM's probability model treats the group order #E(F_p) as a random
integer of size p / exp(delta): delta is the expected log of the
extra smooth part from small primes compared with a random integer,

    delta = sum over small primes l of (E[v_l(order)] - 1/(l-1)) ln l

(GMP-ECM uses delta = 3.134 for Suyama's curves, ECM_EXTRA_SMOOTHNESS in
rho.c, from Kruppa's group-order table).  This script measures E[v_l] of the
group order #E(F_p) for each family on random primes and prints
delta per family and the difference to Suyama's curves, which anchors the
other families to GMP-ECM's 3.134.

  python fam_delta.py bits nsamples [seed]

Families: suyama (-avxecm_curve 0, GMP-ECM param 0), iv (-avxecm_curve 1/2,
GKL family (iv)), j1 (GKL generic Z/2 x Z/4).
"""
import math
import multiprocessing as mp
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import edw_ref as er                         # noqa: E402

bo = er.bo
FAMS = ('suyama', 'iv', 'j1')
SMALL = [q for q in range(2, 200) if bo.is_prime(q)]


def _one(args):
    p, ms = args
    out = {}
    for fam, m in ms.items():
        try:
            A, g, P = er.curve_point(fam, m, p)
        except er.Find:
            out[fam] = None
            continue
        # the group order, not the starting point's order: the index
        # #E / ord(P0) is made of small primes that stage 1 covers anyway,
        # so the rough part (what has to be <= B1, B2) is #E's
        n = er.curve_card(A, g, p)
        out[fam] = [er.v_q(n, q) for q in SMALL]
    return out


def main():
    bits = int(sys.argv[1])
    n = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed)
    jobs = []
    for _ in range(n):
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        while not bo.is_prime(p):
            p += 2
        ms = {'suyama': rng.randrange(6, 1 << 31), 'iv': rng.randrange(1, 1 << 31),
              'j1': rng.randrange(5, 1 << 31)}
        jobs.append((p, ms))
    with mp.Pool() as pool:
        res = pool.map(_one, jobs, chunksize=64)

    print('bits %d samples %d seed %d' % (bits, n, seed))
    delta = {}
    for fam in FAMS:
        rows = [r[fam] for r in res if r[fam] is not None]
        k = len(rows)
        ev = [sum(row[i] for row in rows) / k for i in range(len(SMALL))]
        # per-sample delta, for the standard error
        per = [sum((row[i] - 1.0 / (q - 1)) * math.log(q) for i, q in enumerate(SMALL))
               for row in rows]
        mu = sum(per) / k
        sd = math.sqrt(sum((x - mu) ** 2 for x in per) / (k - 1))
        delta[fam] = (mu, sd / math.sqrt(k), per)
        print('%-7s n=%d  E[v2]=%.4f E[v3]=%.4f E[v5]=%.4f E[v7]=%.4f  delta=%.4f +- %.4f'
              % (fam, k, ev[0], ev[1], ev[2], ev[3], mu, sd / math.sqrt(k)))
    # paired differences against suyama (same primes)
    for fam in FAMS[1:]:
        d = []
        for r in res:
            if r['suyama'] is None or r[fam] is None:
                continue
            a = sum((r[fam][i] - r['suyama'][i]) * math.log(q) for i, q in enumerate(SMALL))
            d.append(a)
        mu = sum(d) / len(d)
        sd = math.sqrt(sum((x - mu) ** 2 for x in d) / (len(d) - 1))
        print('%-7s - suyama: %+.4f +- %.4f  ->  anchored delta %.4f'
              % (fam, mu, sd / math.sqrt(len(d)), 3.134 + mu))


if __name__ == '__main__':
    main()
