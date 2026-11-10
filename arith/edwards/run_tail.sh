#!/bin/bash
# build uedw_tail and run it over the sizes in SPECS ("bits B1" pairs), v2 params
cd /mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8/arith/edwards
gcc -O2 -march=native -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F -DIFMA -DUEDW_NO_PM1 -I../../factor/gmp-ecm \
    uedw_tail.c ../../factor/gmp-ecm/microecm_edw.c -o uedw_tail -lm || exit 1
cd ../..
SPECS=${SPECS:-"36 20 40 27 44 47 48 70 52 85 54 125 56 125 58 165"}
set -- $SPECS
while [ $# -ge 2 ]; do
    arith/edwards/uedw_tail semiprimes_$1bit.dat $2 ${P:-2} 400 | grep -E "B1=|>32|>42|>64|success"
    shift 2
done
