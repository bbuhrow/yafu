#!/usr/bin/env python3
"""Show the prime-power content of microecm's per-level stage-1 scalar s
(as replicated by tiny_chains.prac_level), to define s for new B1 values."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tiny_chains as tc


def primes_upto(n):
    return [q for q in range(2, n + 1) if all(q % d for d in range(2, int(q ** 0.5) + 1))]


def main():
    s70 = tc.parse_steps('microecm.c', 'static void uecm_uprac70(')
    s85 = tc.parse_steps('microecm.c', 'static void uecm_uprac85(')
    for B1 in tc.LEVELS:
        s, adds, dups = tc.prac_level(B1, s70, s85)
        exps = []
        t = s
        for q in primes_upto(4 * B1):
            e = 0
            while t % q == 0:
                t //= q
                e += 1
            if e:
                exps.append('%d^%d' % (q, e) if e > 1 else str(q))
        print('B1 %3d: %3d bits, cofactor %d: %s' % (B1, s.bit_length(), t, ' '.join(exps)))


if __name__ == '__main__':
    main()
