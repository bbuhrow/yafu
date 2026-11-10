#!/bin/bash
# after the P-1 pre-pass change: rebuild yafu, run the harness, the scalar
# timings and the standalone builds that compile microecm_edw.c
set -e
R=/mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8
cd $R
make yafu -j8 2>&1 | grep -E "microecm_edw.c:.*(warning|error)" || true
echo "== harness"
echo 'special()' | ./yafu 2>&1 | grep -E "got"
echo "== scalar timings (v2 column includes the P-1 pre-pass at > 54 bits)"
bash arith/tinyrho/ecm_cmp.sh 4000 | sed -n '1,19p'
echo "== standalone builds"
bash arith/edwards/run_tail.sh > /dev/null && echo run_tail ok
(cd arith/edwards && bash x4_build.sh /tmp/x4chk > /dev/null && echo x4_build ok)
G=$(realpath ../gmp-install/6.2.1)
(cd arith/edwards && gcc -O2 -march=native -I../../factor/gmp-ecm -I../../include -I../../ytools \
    -I$G/include uedw_fails.c ../../factor/gmp-ecm/microecm_edw.c \
    ../../factor/gmp-ecm/micropm1.c $G/lib/libgmp.a -o uedw_fails && echo uedw_fails ok)
arith/edwards/uedw_fails 2 | tail -1
