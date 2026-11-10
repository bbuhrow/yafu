#!/bin/bash
# build and run pm1_eval (micropm1 vs microecm_edw) or pm1_ab (the uecm_edw
# P-1 pre-pass on / off): pm1_build.sh [eval|ab] [args...]
set -e
cd /mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8/arith/tinyrho
R=$(realpath ../..)
G=$(realpath $R/../gmp-install/6.2.1)
D="-O2 -march=native -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F -DIFMA"
INC="-I$R/factor/gmp-ecm -I$R/include -I$R/ytools -I$G/include"
prog=pm1_${1:-eval}; shift || true
gcc $D $INC -o $prog $prog.c $R/factor/gmp-ecm/micropm1.c $R/factor/gmp-ecm/microecm_edw.c $G/lib/libgmp.a -lm
cd $R
arith/tinyrho/$prog "$@"
