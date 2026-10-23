#!/usr/bin/env python3
"""Tune the parameter tables of factor/kmpqs (kmpqs1_param / kmpqs3_param).

usage: tools/kmpqs_tune.py [--variant 1|3] [--rows 5 6 ...] [options]

For each parameter row:
  1. classify a pool of random semiprimes by the row kmpqs uses for them
     (kmpqs_bench -rows; the row depends on the multiplier, not just the size)
  2. coordinate search from the current values: sieve length, FB size (with
     trial division over the whole FB), sieve threshold, number of A factors
     and the size of their pool.  A candidate must factor every input; it is
     kept when it is faster by more than --gain.
  3. re-time the old and new rows on a separate hold-out set and keep the
     new one only if it wins there too.

Prints the tuned tables as C, ready to paste over mpqs_param / mpqs3_param,
and writes them (with all timings) to --out as JSON.

Needs kmpqs_bench (make kmpqs-bench).  Timings are best-of --reps passes, so
run it on an otherwise idle machine.
"""
import argparse, json, math, os, random, subprocess, sys, tempfile

ROWS = {1: 14, 3: 15}
# input sizes whose rows cover each table (row r is roughly bits 4r+41..4r+44
# for kmpqs1 and 4r+89..4r+92 for kmpqs3, shifted by the multiplier)
BITS = {1: range(38, 97), 3: range(86, 149)}
SLMAX = {1: 16384, 3: 16384}
MAXFB = {1: 511, 3: 1023}
MAXADIV = {1: 10, 3: 12}


def run(bench, args):
    r = subprocess.run([bench] + args, capture_output=True, text=True)
    return r.stdout


def current_params(bench):
    tab = {1: {}, 3: {}}
    for line in run(bench, ['-showparams']).split():
        v = [int(x) for x in line.split(',')]
        tab[v[0]][v[1]] = v[2:]
    return tab


def pool(bench, variant, per_bits, seed, tmp):
    """{row: [N, ...]} for random semiprimes over the variant's size range"""
    rows = {}
    for b in BITS[variant]:
        out = run(bench, ['-%d' % variant, '-gen', str(b), str(per_bits),
                          '-seed', str(seed * 1000 + b), '-rows'])
        for line in out.splitlines():
            r, n = line.split()
            rows.setdefault(int(r), []).append(n)
    return rows


def timeit(bench, variant, row, params, path, reps):
    """(seconds, failures) for the inputs in path with row replaced"""
    p = '%d,%d,%s' % (variant, row, ','.join(str(x) for x in params))
    out = run(bench, ['-%d' % variant, '-maxbits', '1000', '-param', p,
                      '-reps', str(reps), path])
    t = fail = None
    for line in out.splitlines():
        if line.startswith('success:'):
            fail = int(line.split()[3])
        elif line.startswith('time:'):
            t = float(line.split()[1])
        elif line.startswith('wrong:') and int(line.split()[1]) > 0:
            fail = 10 ** 6
    if t is None or fail is None:      # rejected or crashed
        return math.inf, 10 ** 6
    return t, fail


def candidates(variant, p):
    """neighbouring parameter rows of p = [nFB, sb, nAdiv, nAdivT, accept, td, sl]"""
    nfb, sb, na, nat, acc, td, sl = p
    out = []
    for f in (0.75, 0.875, 1.125, 1.25, 1.5):
        s = int(round(sl * f / 256.0)) * 256
        if 256 <= s <= SLMAX[variant] and s != sl:
            out.append([nfb, sb, na, nat, acc, td, s])
    for f in (0.8, 0.9, 1.1, 1.25):
        n = int(round(nfb * f))
        if 20 <= n <= MAXFB[variant] and n != nfb:
            out.append([n, sb, na, nat, acc, max(td * n // nfb, 2), sl])
    for d in (-2, -1, 1, 2):
        if 4 <= acc + d <= 40:
            out.append([nfb, sb, na, nat, acc + d, td, sl])
    for d in (-1, 1):
        if 1 <= na + d < nat:
            out.append([nfb, sb, na + d, nat, acc, td, sl])
        if na < nat + d <= MAXADIV[variant]:
            out.append([nfb, sb, na, nat + d, acc, td, sl])
    return out


def tune_row(bench, variant, row, p0, train, hold, reps, gain, passes, log):
    best, (bt, bf) = list(p0), timeit(bench, variant, row, p0, train, reps)
    log('  start %s: %.4f s, %d failures' % (best, bt, bf))
    base_fail = bf
    for ps in range(passes):
        improved = False
        for c in candidates(variant, best):
            # paired: time the candidate and the current best back to back,
            # twice, so drift in machine speed cancels
            t, f = timeit(bench, variant, row, c, train, reps)
            if f > base_fail:
                continue
            bt, _ = timeit(bench, variant, row, best, train, reps)
            if t >= bt * (1 - gain):
                continue
            t2, f2 = timeit(bench, variant, row, c, train, reps)
            bt2, _ = timeit(bench, variant, row, best, train, reps)
            if f2 <= base_fail and t2 < bt2 * (1 - gain / 2):
                log('    %s: %.4f s vs %.4f s  (%.1f%% better)'
                    % (c, t2, bt2, 100 * (1 - t2 / bt2)))
                best, bt, improved = c, t2, True
        if not improved:
            break
    if best == list(p0):
        return best, None
    t_old, f_old = timeit(bench, variant, row, p0, hold, reps)
    t_new, f_new = timeit(bench, variant, row, best, hold, reps)
    log('  hold-out: old %.4f s (%d fail), new %.4f s (%d fail)' % (t_old, f_old, t_new, f_new))
    if f_new <= f_old and t_new < t_old * (1 - gain):
        return best, (t_old, t_new)
    return list(p0), None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bench', default='./kmpqs_bench')
    ap.add_argument('--variant', type=int, nargs='+', default=[1, 3])
    ap.add_argument('--rows', type=int, nargs='+', help='rows to tune (default all)')
    ap.add_argument('--count', type=int, default=1000, help='training inputs per row')
    ap.add_argument('--holdout', type=int, default=1000, help='hold-out inputs per row')
    ap.add_argument('--reps', type=int, default=3)
    ap.add_argument('--budget', type=float, default=0.5,
                    help='target seconds for one pass over the training inputs')
    ap.add_argument('--gain', type=float, default=0.015,
                    help='minimum relative speedup to accept a change')
    ap.add_argument('--passes', type=int, default=3)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default='kmpqs_tune.json')
    a = ap.parse_args()

    logf = open(os.path.splitext(a.out)[0] + '.log', 'a')

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + '\n')
        logf.flush()

    tab = current_params(a.bench)
    result = {'tables': {}, 'timings': {}}
    tmp = tempfile.mkdtemp(prefix='kmpqs_tune_')
    for v in a.variant:
        log('=== kmpqs%d: classifying inputs' % v)
        # each row spans about 4 input sizes; ask for a margin
        per_bits = max(10, (a.count + a.holdout) // 2)
        P = pool(a.bench, v, per_bits, a.seed, tmp)
        rows = a.rows if a.rows else range(ROWS[v])
        for r in rows:
            ns = P.get(r, [])
            random.Random(r).shuffle(ns)
            if len(ns) < 20:
                log('row %d: only %d inputs, skipped' % (r, len(ns)))
                continue
            k = min(a.count, len(ns) // 2)
            train = os.path.join(tmp, 'k%d_r%d_train.txt' % (v, r))
            hold = os.path.join(tmp, 'k%d_r%d_hold.txt' % (v, r))
            # keep one timing pass near a.budget seconds
            open(train, 'w').write('\n'.join(ns[:k]) + '\n')
            t0, _ = timeit(a.bench, v, r, tab[v][r], train, 1)
            if t0 > a.budget:
                k = max(60, int(k * a.budget / t0))
                open(train, 'w').write('\n'.join(ns[:k]) + '\n')
            open(hold, 'w').write('\n'.join(ns[k:k + min(a.holdout, k)]) + '\n')
            log("row %d (%d train, %d hold-out inputs)" % (r, k, min(a.holdout, k, len(ns) - k)))
            new, tim = tune_row(a.bench, v, r, tab[v][r], train, hold, a.reps,
                                a.gain, a.passes, log)
            if tim:
                log('  row %d: %s -> %s  (hold-out %.1f%% faster)'
                    % (r, tab[v][r], new, 100 * (1 - tim[1] / tim[0])))
                result['timings']['%d/%d' % (v, r)] = tim
            else:
                log('  row %d: unchanged' % r)
            tab[v][r] = new
        result['tables'][v] = [tab[v][r] for r in range(ROWS[v])]
        log('\nstatic uint16_t mpqs%s_param[%d][7]={' % ('' if v == 1 else '3', ROWS[v]))
        for r in range(ROWS[v]):
            log(' { %s }%s' % (', '.join(str(x) for x in tab[v][r]),
                               ',' if r + 1 < ROWS[v] else ''))
        log('};')
        json.dump(result, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
