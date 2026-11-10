#!/bin/bash
# Build the end-to-end benchmark (tiny_e2e.c, with the Edwards and the AVX2
# x4 rows) and the x4 test drivers, per ISA, from the repo sources.
#   x4_build.sh OUTDIR
# avx2: -mavx2 only (FMA only for the x4 files, as in the Makefile);
# ifma / f: AVX-512 builds (x8 rows), which also run the x4 code.
# X4FLAGS (env): extra flags for the x4 files, e.g. -DUEDW4_W=1.
set -e
OUT=$(realpath -m "$1")
mkdir -p "$OUT"
R=$(realpath ../..)
S=$R/factor/gmp-ecm
G=$(realpath $R/../gmp-install/6.2.1)
INC="-I$S -I$R/ytools -I$R/include -I$G/include"
for isa in avx2 ifma f; do
    case $isa in
        avx2) D="-O2 -mavx2 -mbmi -mbmi2 -DUSE_AVX2 -DUSE_BMI2";;
        ifma) D="-O2 -march=native -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F -DIFMA";;
        f)    D="-O2 -march=native -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F";;
    esac
    O=$OUT/o_$isa; mkdir -p $O
    for src in microecm micropm1 tinyecm microecm_edw tinyecm_edw; do
        gcc $D $INC -c $S/$src.c -o $O/$src.o
    done
    gcc $D -mfma $X4FLAGS $INC -c $S/microecm_edw_x4.c -o $O/microecm_edw_x4.o
    if [ -f $S/tinyecm_edw_x4.c ]; then
        gcc $D -mfma $X4FLAGS $INC -c $S/tinyecm_edw_x4.c -o $O/tinyecm_edw_x4.o
    fi
    gcc $D $INC -c $R/arith/monty.c -o $O/monty.o
    gcc $D $INC -DWITH_EDW -DWITH_X4 -o $OUT/e2e_$isa tiny_e2e.c $O/*.o $G/lib/libgmp.a -lm
done
echo built in $OUT
