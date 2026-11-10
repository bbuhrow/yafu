#!/bin/bash
# microecm_edw timings on the same sizes, for comparison with tinyrho:
# scalar (mode 0) and the x8 list (mode 1, n <= 52 bits)
set -e
cd /mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8/arith/edwards
G=$(realpath ../../../gmp-install/6.2.1)
gcc -O2 -march=native -mfma -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F -DIFMA -I../../factor/gmp-ecm -I../../include -I../../ytools -I$G/include \
    uedw_v2bench.c ../../factor/gmp-ecm/microecm_edw.c ../../factor/gmp-ecm/microecm_edw_x4.c \
    ../../factor/gmp-ecm/micropm1.c \
    $G/lib/libgmp.a -lm -o /tmp/v2b
echo "== scalar"; /tmp/v2b ${1:-4000} 2 1 0 32 64 0
echo "== x8 list"; /tmp/v2b ${1:-4000} 2 1 0 32 52 1
