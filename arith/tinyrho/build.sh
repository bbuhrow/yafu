#!/bin/bash
# Build rho_test per ISA: ifma (AVX-512 IFMA), f (AVX-512F, no IFMA),
# avx2 (no AVX-512: the list routines fall back to the scalar code).
#   build.sh [OUTDIR]       (default: arith/tinyrho)
set -e
cd "$(dirname "$0")"
OUT=$(realpath -m "${1:-.}")
R=$(realpath ../..)
G=$(realpath $R/../gmp-install/6.2.1)
INC="-I$R/factor/gmp-ecm -I$R/arith/limb -I$G/include"
for isa in ifma f avx2; do
    case $isa in
        ifma) D="-O2 -march=native -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F -DIFMA";;
        f)    D="-O2 -march=native -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F";;
        avx2) D="-O2 -mavx2 -mbmi -mbmi2 -DUSE_AVX2 -DUSE_BMI2";;
    esac
    gcc $D -Wall $INC -o $OUT/rho_test_$isa rho_test.c $R/factor/gmp-ecm/tinyrho.c $G/lib/libgmp.a -lm
done
echo built in $OUT
