#!/usr/bin/env python3
"""Compare b1b2_study's model with uedw_lvlbench measurements (single
factor, fixed level): b1b2_check.py b1b2_lvlbench.txt"""
import re
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b1b2_study as st         # noqa: E402
import b1b2_emit as em          # noqa: E402

V2 = {B1: (pol, m) for B1, pol, m in em.LEVELS2}
rows = {}
print('  N set  B1 |  real ns  model ns  ratio | real curves  model curves')
for line in open(sys.argv[1]):
    m = re.match(r'N (\d+) set (\d) B1\s+(\d+):\s+(\d+) ns per find\s+(\d+) ns per input\s+curves '
                 r'([\d.]+)\s+found ([\d.]+)', line)
    if not m:
        continue
    N, sset, B1 = int(m.group(1)), int(m.group(2)), int(m.group(3))
    real_t, real_c = float(m.group(5)), float(m.group(6))
    if N not in rows:
        rows[N] = st.load(N)
    if sset == 1:
        pol, mult = 'lvl', st.legacy_mult(B1)
    else:
        pol, mult = V2[B1]
    s = st.get_scalar(B1, pol)
    t, fr, nc = st.eval_one(rows[N], s, st.stage2(B1, mult, pol), st.t1_cost(B1, pol), 'one')
    # per input (the bench reports per input; the model per find, ~equal here)
    print('%3d %3d %3d | %8.0f %8.0f %6.3f | %6.2f %6.2f' % (N, sset, B1, real_t, t * fr,
                                                          real_t / (t * fr), real_c, nc))
