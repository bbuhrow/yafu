#!/usr/bin/env python3
"""Test-sieve one GNFS parameter set, bayestune style (factor/tune.c bt_nfs_ref):
start-up measured with 1 q per thread, a 100-q pilot at q0 sizes the blocks
(~BLOCK_SEC per thread) and the span (1.2 * need / y0), then one block per
thread at q spread over the span, all threads at once; sieve time integrated
linearly in q (bt_nfs_integrate) until 'need' raw relations.

usage: testsieve_params.py POLY TAG rlim alim lpbr lpba mfbr mfba rlam alam I q0 need [threads] [blocksec]
"""
import os, subprocess, sys, time, shutil

SIEVEDIR = os.environ.get("SIEVEDIR", "factor/lasieve5_64/bin/avx512")
RES = os.path.expanduser("~/ts_work/results2.txt")

def spread(siever, wd, job, qv, rng):
    procs = []
    for i, q in enumerate(qv):
        out = os.path.join(wd, "out.%d" % i)
        if os.path.exists(out):
            os.remove(out)
        t0 = time.time()
        p = subprocess.Popen([siever, "-f", str(int(q)), "-c", str(int(rng)), "-o", out,
                              "-n", str(i), "-a", job], cwd=wd,
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        procs.append((p, t0, out))
    secs, rels = [], []
    for p, t0, out in procs:
        p.wait()
        secs.append(time.time() - t0)
        n = 0
        if os.path.exists(out):
            with open(out) as f:
                n = sum(1 for l in f if l and l[0] != "#")
            os.remove(out)
        rels.append(n)
    return secs, rels

def integrate(qs, y, s, q0, need):
    q, rels, t = q0, 0.0, 0.0
    dq = max(1.0, need / y[0] / 4000.0)
    it = 0
    while rels < need and it < 1000000:
        it += 1
        qm = q + dq / 2
        if qm <= qs[0]:
            yy, ss = y[0], s[0]
        elif qm >= qs[-1]:
            yy, ss = y[-1], s[-1]
        else:
            k = 1
            while qm >= qs[k]:
                k += 1
            f = (qm - qs[k - 1]) / (qs[k] - qs[k - 1])
            yy = y[k - 1] + (y[k] - y[k - 1]) * f
            ss = s[k - 1] + (s[k] - s[k - 1]) * f
        if yy <= 0:
            break
        rels += yy * dq
        t += ss * dq
        q += dq
    return t, q

def main():
    a = sys.argv[1:]
    poly, tag = a[0], a[1]
    rlim, alim, lpbr, lpba, mfbr, mfba = [int(float(x)) for x in a[2:8]]
    rlam, alam = float(a[8]), float(a[9])
    I = int(a[10])
    q0 = float(a[11])
    need = float(a[12])
    T = int(a[13]) if len(a) > 13 else 6
    blocksec = float(a[14]) if len(a) > 14 else 120.0
    siever = os.path.abspath(os.path.join(SIEVEDIR, "gnfs-lasieve4I%de" % I))

    wd = os.path.expanduser("~/ts_work/" + tag)
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    keys = ("n", "skew", "c0", "c1", "c2", "c3", "c4", "c5", "c6", "Y0", "Y1")
    lines = [l for l in open(poly) if l.split(":")[0].strip() in keys]
    job = os.path.join(wd, "job")
    with open(job, "w") as f:
        f.writelines(lines)
        f.write("rlim: %d\nalim: %d\nlpbr: %d\nlpba: %d\nmfbr: %d\nmfba: %d\n"
                "rlambda: %.3f\nalambda: %.3f\n" % (rlim, alim, lpbr, lpba, mfbr, mfba, rlam, alam))

    secs, _ = spread(siever, wd, job, [q0 + i for i in range(T)], 1)
    startup = sum(secs) / T
    secs, rels = spread(siever, wd, job, [q0 + 100 * i for i in range(T)], 100)
    secq = sum(max(x - startup, 0.2 * x) for x in secs) / (100.0 * T)
    y0 = sum(rels) / (100.0 * T)
    R = max(50, min(200000, 10 * round(blocksec / secq / 10)))
    S = 1.2 * need / y0
    qs = [q0 + S * i / (T - 1) for i in range(T)]
    secs, rels = spread(siever, wd, job, qs, R)
    y, s, detail = [], [], []
    for i in range(T):
        busy = max(secs[i] - startup, 0.5 * secs[i])
        y.append(rels[i] / R)
        s.append(busy / R / T)
        detail.append("%.1fM:%.2f/%.0f" % (qs[i] / 1e6, rels[i] / R, rels[i] / busy))
    sieve, qend = integrate(qs, y, s, q0, need)
    line = ("%s | lim %.2fM/%.2fM lpb %d/%d mfb %d/%d lam %.3f/%.3f %de q0 %.2fM need %.1fM"
            " | sieve %.0f s (%.1f h) q %.1fM-%.1fM | startup %.1f R %d | %s" % (
            tag, rlim / 1e6, alim / 1e6, lpbr, lpba, mfbr, mfba, rlam, alam, I, q0 / 1e6,
            need / 1e6, sieve, sieve / 3600, q0 / 1e6, qend / 1e6, startup, R, " ".join(detail)))
    print(line, flush=True)
    with open(RES, "a") as f:
        f.write(line + "\n")
    shutil.rmtree(wd, ignore_errors=True)

main()
