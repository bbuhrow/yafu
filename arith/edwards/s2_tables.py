#!/usr/bin/env python3
"""Stage-2 pairing tables for microecm's D = 60 x-only stage 2.

The stage 2 walks giant points A = A0, A0 + 120, ... (A = 30 mod 60, or the
special starts of the existing levels) and for each A forms one cross term
per baby b in {1, 7, 11, ..., 59} (coprime to 30), which catches a prime r
when r = A - b or r = A + b.  table(B1, B2) builds the natural table: per A,
every b with A - b or A + b a prime in (B1, B2].

  s2_tables.py cmp      compare with microecm's tables (covered primes, size)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BABY = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59]


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def table(B1, B2, A0=None):
    """[(A, [b, ...]), ...]: giant steps from A0 (default: the first
    A = 30 mod 120 window reaching above B1) covering the primes in
    (B1, B2] (primes of 2, 3, 5 excluded: never in the windows)"""
    P = set(q for q in primes_upto(B2 + 60) if B1 < q <= B2 and q > 5)
    if A0 is None:
        A0 = 30
        while A0 + 59 <= B1:
            A0 += 120
    out = []
    A = A0
    while A - 59 <= B2:
        bs = [b for b in BABY if (A - b) in P or (A + b) in P]
        out.append((A, bs))
        A += 120
    return out


def covered(tab):
    c = set()
    for A, bs in tab:
        for b in bs:
            c.add(A - b)
            c.add(A + b)
    return c


def cost(tab):
    """(cross terms, giant steps)"""
    return sum(len(bs) for _, bs in tab), len(tab)


def main():
    import tiny_curves as tcv
    for B1 in (27, 47, 70, 85, 125, 165, 205):
        pr = tcv.stage2_pairs(B1)
        cov = set()
        for A, b in pr:
            cov.add(A - b)
            cov.add(A + b)
        As = sorted(set(A for A, _ in pr))
        Bmax = max(cov)
        mult = 50 if B1 >= 125 else 25
        B2 = mult * B1
        P = [q for q in primes_upto(B2 + 200) if B1 < q]
        miss = [q for q in P if q <= B2 and q not in cov]
        extra = [q for q in P if q > B2 and q in cov]
        nat = table(B1, B2)
        nc, ng = cost(nat)
        print('B1 %3d B2 %5d (x%d): microecm %4d crosses %3d giants, A %d..%d, covers to %d, '
              'misses %d primes <= B2 (first %s), %d beyond | natural %4d crosses %3d giants'
              % (B1, B2, mult, len(pr), len(As), As[0], As[-1], Bmax, len(miss), miss[:4],
                 len(extra), nc, ng))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'cmp':
        main()
