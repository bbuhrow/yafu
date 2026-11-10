"""
Phase 0 stage-1 cost model: yafu's PRAC on Montgomery curves vs a = -1
twisted Edwards with signed windows (arith/avxecm_edwards_plan.md).

  python edw_model.py [B1 ...]

Operation counts are exact for yafu's stage-1 scalar (prime powers q^k < B1):
  * PRAC: replays vec_prac (avxecm.c, ORIG_PRAC off: conditions 3, 4, 5, 9;
    v chosen by lucas_cost with ADD = 5.5, DUP = 4.5) and counts its steps;
  * Edwards: width-w NAF of the whole scalar (one batch; the C driver chunks
    per prime range, which changes nothing at these lengths), table of
    2^(w-2) cached odd multiples.
Each step is priced with ops_<isa>.txt (factor/avx-ecm/bench: ae_bench -e,
ns per 8-lane call of the composed operations, as the curve code calls the
kernels):
  PRAC condition 3 / final add   mont_add
  PRAC conditions 4, 5, 9        mont_add + mont_dup - add - sub
                                 (the doubling reuses the add's sum/diff)
  PRAC start (per prime power)   mont_dup + 6 copies (~ 6 add)
  doublings for 2^k              mont_dup
  Edwards                        edw_dbl / edw_dblt / edw_add per digit;
                                 table: edw_dblt + (2^(w-2) - 1)(edw_add + mul)
                                 + 2^(w-2) (mul + 3 add)  (T and cached form)
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COLS = ['mul', 'sqr', 'add', 'sub', 'mont_add', 'mont_dup', 'edw_dbl', 'edw_dblt', 'edw_add']
VAL = [0.61803398874989485, 0.72360679774997897, 0.58017872829546410,
       0.63283980608870629, 0.61242994950949500, 0.62018198080741576,
       0.61721461653440386, 0.61834711965622806, 0.61791440652881789,
       0.61807966846989581]
ADD, DUP = 5.5, 4.5


def primes_below(n):
    sieve = bytearray([1]) * n
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, n, i)))
    return [i for i in range(n) if sieve[i]]


def lucas_cost(n, v):
    d = n
    r = int(d * v + 0.5)
    if r >= n:
        return ADD * n
    d = n - r
    e = 2 * r - n
    c = DUP + ADD
    while d != e:
        if d < e:
            d, e = e, d
        if (d + 3) // 4 <= e:
            d -= e
            c += ADD
        elif (d + e) % 2 == 0:
            d = (d - e) // 2
            c += ADD + DUP
        elif d % 2 == 0:
            d //= 2
            c += ADD + DUP
        else:
            e //= 2
            c += ADD + DUP
    return c


def prac_counts(q):
    """(cond-3 adds, add+dup pairs) of vec_prac(q); plus 1 start dup and
    1 final add, counted by the caller"""
    best, cmin = 0, 1e300
    for i, v in enumerate(VAL):
        cst = lucas_cost(q, v)
        if cst < cmin:
            cmin, best = cst, i
    r = int(q * VAL[best] + 0.5)
    d = q - r
    e = 2 * r - q
    n3 = npair = 0
    while d != e:
        if d < e:
            d, e = e, d
        if (d + 3) // 4 <= e:
            d -= e
            n3 += 1
        elif (d + e) % 2 == 0:
            d = (d - e) // 2
            npair += 1
        elif d % 2 == 0:
            d //= 2
            npair += 1
        else:
            e //= 2
            npair += 1
    return n3, npair


def stage1_powers(B1):
    """yafu's stage-1 prime powers: 2^k < B1 as doublings, q^j with q^j < B1
    as repeated vec_prac(q)"""
    k2 = 0
    q = 2
    while q < B1:
        k2 += 1
        q *= 2
    pw = []
    for q in primes_below(B1)[1:]:
        c = 1
        while True:
            pw.append(q)
            c *= q
            if c * q >= B1:
                break
    return k2, pw


def scalar(B1):
    k2, pw = stage1_powers(B1)
    # product tree for speed
    vals = pw[:]
    while len(vals) > 1:
        vals = [vals[i] * vals[i + 1] if i + 1 < len(vals) else vals[i]
                for i in range(0, len(vals), 2)]
    return vals[0] << k2, k2, pw


def wnaf_counts(k, w):
    """(dbl, dbl_t, add) for the main loop of a width-w NAF of k"""
    half, full = 1 << (w - 1), 1 << w
    nd = nz = 0
    while k:
        if k & 1:
            r = k & (full - 1)
            if r >= half:
                r -= full
            k -= r
            nz += 1
        nd += 1
        k >>= 1
    nd -= 1                 # top digit: no doubling before it
    nz -= 1                 # the top digit initializes R (no add)
    return nd - nz, nz, nz


def load_ops(isa):
    rows = {}
    with open(os.path.join(HERE, 'ops_%s.txt' % isa)) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 11 and parts[0].isdigit():
                rows[int(parts[0])] = dict(zip(COLS, map(float, parts[2:])))
    return rows


def main(B1s):
    isas = [i for i in ('ifma', 'noifma', 'avx2') if os.path.exists(os.path.join(HERE, 'ops_%s.txt' % i))]
    ops = {i: load_ops(i) for i in isas}
    for B1 in B1s:
        s, k2, pw = scalar(B1)
        bits = s.bit_length()
        n3 = npair = 0
        cache = {}
        for q in pw:
            if q not in cache:
                cache[q] = prac_counts(q)
            a, b = cache[q]
            n3 += a
            npair += b
        nstart = nfinal = len(pw)
        print('B1 = %d: scalar %d bits; PRAC: %d prime powers, %.3f adds + %.3f doublings per bit'
              % (B1, bits, len(pw), (n3 + npair + nfinal) / bits, (npair + nstart + k2) / bits))
        wc = {}
        for w in range(4, 14):
            wc[w] = wnaf_counts(s, w)
            nd, ndt, na = wc[w]
            print('  wNAF w=%2d: %.4f dbl + %.4f dbl_t + %.4f add per bit (density %.4f), table %d points'
                  % (w, nd / bits, ndt / bits, na / bits, na / bits, 1 << (w - 2)))
        for isa in isas:
            print('  %s: ns per bit of s (8 curves); ratio = PRAC / Edwards' % isa)
            print('     S  bits    PRAC   Edw w=6   w=8  w=10  w=12 | best w  ratio  table MB(best w)')
            for S, o in sorted(ops[isa].items()):
                pair = o['mont_add'] + o['mont_dup'] - o['add'] - o['sub']
                tprac = (n3 * o['mont_add'] + npair * pair + nstart * (o['mont_dup'] + 6 * o['add'])
                         + nfinal * o['mont_add'] + k2 * o['mont_dup']) / bits
                ted = {}
                for w, (nd, ndt, na) in wc.items():
                    nt = 1 << (w - 2)
                    tab = o['edw_dblt'] + (nt - 1) * (o['edw_add'] + o['mul']) + nt * (o['mul'] + 3 * o['add'])
                    ted[w] = (nd * o['edw_dbl'] + ndt * o['edw_dblt'] + na * o['edw_add'] + tab) / bits
                wb = min(ted, key=ted.get)
                mb = (1 << (wb - 2)) * 4 * S * 64 / 2 ** 20
                print('    %2d %5d %7.0f %7.0f %5.0f %5.0f %5.0f |  %2d  %5.3f  %6.2f'
                      % (S, 52 * S, tprac, ted[6], ted[8], ted[10], ted[12], wb, tprac / ted[wb], mb))


if __name__ == '__main__':
    B1s = [int(float(a)) for a in sys.argv[1:]] or [100000, 1000000]
    main(B1s)
