#!/bin/bash
# Check that every kernel ISA of factor/kmpqs gives identical results, and
# time each one.
#
# usage: tools/kmpqs_isa_check.sh [sizes...]      (default 64 80 96 100 120)
#
# Builds kmpqs_bench_stat for each ISA (avx512 avx2 sse41 sse2 scalar, the
# ones the build machine can run) into $WORK, generates semiprimes per size
# (fixed seed), and compares the "Stat:" counters of every ISA against the
# scalar build.  The kernels compute exactly what lasieve's C code does, so
# the counters must be identical.
#
# Environment:
#   WORK        work directory           (default ~/kmpqs_isa)
#   GMP_PREFIX  GMP install prefix       (default <src>/../gmp-install/6.2.1)
#   COUNT       inputs per size          (default 1000)
#   REPS        timing passes, best one reported (default 3)
#   CC          compiler                 (default gcc)

set -e

SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=${WORK:-$HOME/kmpqs_isa}
GMP_PREFIX=${GMP_PREFIX:-$(cd "$SRC/.." && pwd)/gmp-install/6.2.1}
SIZES=${*:-64 80 96 100 120}
CC=${CC:-gcc}
mkdir -p "$WORK"

KSRC="factor/kmpqs/kmpqs.c factor/kmpqs/kmpqs1.c factor/kmpqs/kmpqs3.c factor/kmpqs/kmpqs_bench.c
      top/aprcl/tinyprp.c arith/monty.c ytools/ytools.c"
INC="-I$SRC/include -I$SRC/ytools -I$SRC/top/aprcl -I$GMP_PREFIX/include"

declare -A FLAGS=(
	[avx512]="-march=icelake-client"
	[avx2]="-march=haswell"
	[sse41]="-march=nehalem"
	[sse2]="-march=x86-64"
	[scalar]="-march=x86-64 -DKMPQS_NO_SIMD"
)
ISAS="scalar sse2 sse41 avx2 avx512"
grep -q avx512bw /proc/cpuinfo || ISAS="scalar sse2 sse41 avx2"

for isa in $ISAS; do
	for stat in 0 1; do
		out="$WORK/kb_${isa}"; extra=""
		if [ $stat = 1 ]; then out="${out}_stat"; extra="-DMPQS_STAT -DMPQS3_STAT"; fi
		(cd "$SRC" && $CC -O2 ${FLAGS[$isa]} $extra -w $INC $KSRC \
			"$GMP_PREFIX/lib/libgmp.a" -lm -lpthread -o "$out")
	done
done

cd "$WORK"
status=0
printf "%5s" bits; for isa in $ISAS; do printf " %9s" "$isa"; done; printf "   counters\n"
for b in $SIZES; do
	in="isa_${b}.txt"
	./kb_scalar -gen "$b" "${COUNT:-1000}" -seed $((b + 1)) -dump > "$in"
	if [ "$b" -le 96 ]; then v=-1; else v=-3; fi
	ref=$(./kb_scalar_stat $v "$in" | grep -E '^(success|Stat|  ff|wrong)')
	verdict=same
	printf "%5s" "$b"
	for isa in $ISAS; do
		c=$(./kb_${isa}_stat $v "$in" | grep -E '^(success|Stat|  ff|wrong)')
		[ "$c" = "$ref" ] || { verdict="DIFF($isa)"; status=1; }
		t=$(./kb_${isa} $v -reps ${REPS:-3} "$in" | grep '^time:' | awk '{print $2}')
		printf " %9s" "$t"
	done
	printf "   %s\n" "$verdict"
done
exit $status
