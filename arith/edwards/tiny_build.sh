#!/bin/bash
# Build the micro/tinyecm benchmarks against a given copy of the sources.
#   tiny_build.sh SRCDIR OUTDIR
# SRCDIR holds microecm.c, tinyecm.c and monty.c (e.g. the repo's, or
# copies from an older commit); headers come from the repo.  Builds, per
# ISA: e2e (tiny_e2e.c through the public entry points) and the stage-level
# benches (tiny_bench*.c, which #include the .c files).  Run from arith/edwards.
set -e
SRC=$(realpath "$1"); OUT=$(realpath -m "$2")
mkdir -p "$OUT"
R=$(realpath ../..)
G=$(realpath $R/../gmp-install/6.2.1)
INC="-I$SRC -I$R/factor/gmp-ecm -I$R/arith/limb -I$R/ytools -I$R/include -I$G/include"
# the 64/128-bit kernels live in the limb library now (arith/limb)
LIMB="$R/arith/limb/limb1.c $R/arith/limb/limb2.c"
BASE="-O2 -march=native -DUSE_AVX2 -DUSE_BMI2"
for isa in ifma f; do
    case $isa in
        ifma) D="-DUSE_AVX512F -DIFMA";;
        f)    D="-DUSE_AVX512F";;
    esac
    # stage-level benches include the sources through a shim directory
    S=$OUT/shim_$isa; mkdir -p $S/x/y
    cp tiny_bench.c tiny_bench8.c tiny_bench104.c tiny_bench128.c edw_tiny_chains.h $S/x/y/
    mkdir -p $S/factor/gmp-ecm
    cp $SRC/microecm.c $SRC/tinyecm.c $S/factor/gmp-ecm/
    cp $R/factor/gmp-ecm/*.h $S/factor/gmp-ecm/
    (cd $S/x/y &&
     gcc $BASE $D $INC -o $OUT/bench64_$isa tiny_bench.c -lm &&
     gcc $BASE $D $INC -o $OUT/bench8_$isa tiny_bench8.c -lm &&
     gcc $BASE $D $INC -o $OUT/bench104_$isa tiny_bench104.c $SRC/monty.c $LIMB $G/lib/libgmp.a -lm &&
     gcc $BASE $D $INC -o $OUT/bench128_$isa tiny_bench128.c $SRC/monty.c $LIMB $G/lib/libgmp.a -lm)
    gcc $BASE $D $INC -c $SRC/microecm.c -o $OUT/m_$isa.o
    gcc $BASE $D $INC -c $SRC/tinyecm.c -o $OUT/t_$isa.o
    gcc $BASE $D $INC -c $SRC/monty.c -o $OUT/mo_$isa.o
    gcc $BASE $D $INC -o $OUT/e2e_$isa tiny_e2e.c $OUT/m_$isa.o $OUT/t_$isa.o $OUT/mo_$isa.o \
        $LIMB $G/lib/libgmp.a -lm
done
echo built in $OUT
