#!/usr/bin/env python3
"""Summarize tools/bench_matrix.sh CSV output.

usage: bench_summary.py results.csv [--base LABEL] [--groups "a,b;c,d"]

Prints mean QS seconds per (input, NB, variant) with % vs the base variant
(default: the first variant in the file).  Also checks the work counters:
variants listed in the same --groups group (comma-separated labels; groups
separated by ';') must have identical polys / touched / blocks / full /
partial for each (input, NB, rep); with no --groups, all variants are
expected to match.  Any run that didn't factor is flagged.
"""
import csv
import sys
import argparse
from collections import defaultdict, OrderedDict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--base")
    ap.add_argument("--groups")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.csv)))
    if not rows:
        print("no rows")
        return 1
    variants = list(OrderedDict.fromkeys(r["variant"] for r in rows))
    base = a.base or variants[0]
    if a.groups:
        groups = [g.split(",") for g in a.groups.split(";")]
    else:
        groups = [variants]

    t = defaultdict(list)
    for r in rows:
        if r["qs_sec"]:
            t[(r["input"], r["nb"], r["variant"])].append(float(r["qs_sec"]))

    keys = list(OrderedDict.fromkeys((r["input"], r["nb"]) for r in rows))
    w = max(10, max(len(v) for v in variants) + 2)
    print("%-6s %-4s " % ("input", "NB") + "".join("%*s" % (w + 8, v) for v in variants))
    for (inp, nb) in keys:
        b = t.get((inp, nb, base))
        bm = sum(b) / len(b) if b else None
        cells = []
        for v in variants:
            x = t.get((inp, nb, v))
            if not x:
                cells.append("%*s" % (w + 8, "-"))
                continue
            m = sum(x) / len(x)
            pct = "" if (bm is None or v == base) else " (%+.1f%%)" % (100 * (m - bm) / bm)
            cells.append("%*s" % (w + 8, "%.2f%s" % (m, pct)))
        print("%-6s %-4s " % (inp, nb) + "".join(cells))

    bad = 0
    work = defaultdict(dict)
    for r in rows:
        if r["factored"] != "1":
            print("NOT FACTORED: %s %s NB%s rep%s" % (r["variant"], r["input"], r["nb"], r["rep"]))
            bad += 1
        work[(r["input"], r["nb"], r["rep"])][r["variant"]] = tuple(
            r[k] for k in ("polys", "touched", "blocks", "full", "partial"))
    for g in groups:
        for k, d in work.items():
            vals = {v: d[v] for v in g if v in d}
            if len(set(vals.values())) > 1:
                bad += 1
                print("WORK MISMATCH %s NB%s rep%s:" % k)
                for v, x in vals.items():
                    print("    %-14s polys %s touched %s blocks %s full %s partial %s" % ((v,) + x))
    if bad == 0:
        print("work counters: identical within groups; all runs factored")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
