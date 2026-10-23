#!/bin/bash
# Validate and benchmark factor/kmpqs against the lasieve5_64 originals.
#
# usage: tools/kmpqs_compare.sh [sizes...]
#
# Builds, in $WORK (default ~/kmpqs_ref):
#   lasieve_asm/    factor/lasieve5_64 as shipped (m4/GAS asm, Linux x86-64)
#   lasieve_noasm/  the same with every ASM_MPQS*, HAVE_XMM_MUL and
#                   A64_STYLE_TD switch in asm/mpqs-config.h turned off
#   kmpqs_bench, kmpqs_bench_stat   from this tree (make kmpqs-bench)
#
# For each size it generates semiprimes with kmpqs_bench -gen (fixed seed),
# then
#   - diffs kmpqs_bench_stat's counters against lasieve_noasm's mpqstest /
#     mpqs3test.  These must match exactly up to 100 bits: kmpqs is lasieve's
#     C code.  Above ~104 bits lasieve's mpqs3 corrupts relations (see
#     mpqs3_convert in factor/kmpqs/kmpqs3.c) and mostly fails, so larger
#     sizes are only checked for correctness, not compared.
#   - times lasieve_asm and kmpqs_bench on the same inputs.
#
# Sizes <= 96 use mpqs/kmpqs1, larger sizes mpqs3/kmpqs3.  Default sizes:
#   50 64 80 90 96 100 110 120
#
# Environment:
#   WORK        work directory           (default ~/kmpqs_ref)
#   GMP_PREFIX  GMP install prefix       (default <src>/../gmp-install/6.2.1)
#   COUNT       inputs per size          (default 2000; 500 above 100 bits)
#   KMPQS_ONLY  =1: skip the lasieve builds and comparisons, just time kmpqs
#   REPS        timing runs per program and size, best one reported (default 5)

set -e

SRC=$(cd "$(dirname "$0")/.." && pwd)
WORK=${WORK:-$HOME/kmpqs_ref}
GMP_PREFIX=${GMP_PREFIX:-$(cd "$SRC/.." && pwd)/gmp-install/6.2.1}
SIZES=${*:-50 64 80 90 96 100 110 120}
mkdir -p "$WORK"

# --- kmpqs ------------------------------------------------------------------
make -C "$SRC" kmpqs-bench GMP_PREFIX="$GMP_PREFIX" >/dev/null
cp "$SRC/kmpqs_bench" "$SRC/kmpqs_bench_stat" "$WORK/"

# --- lasieve references -----------------------------------------------------
build_lasieve() {    # $1 = dir name, $2 = 1 to strip the asm switches
	local d="$WORK/$1"
	if [ -x "$d/mpqstest" ] && [ -x "$d/mpqs3test" ]; then return; fi
	rm -rf "$d"
	cp -r "$SRC/factor/lasieve5_64" "$d"
	rm -f "$d"/*.o "$d"/asm/*.o "$d"/asm/*.a
	if [ "$2" = 1 ]; then
		sed -i -E 's@^#define (ASM_MPQS|HAVE_XMM_MUL|A64_STYLE_TD)@// #define \1@' \
			"$d/asm/mpqs-config.h"
	fi
	make -C "$d/asm" liblasieve.a GMP_INCDIR="$GMP_PREFIX/include" \
		DETECTED_OS=Linux >/dev/null 2>&1
	make -C "$d" mpqstest mpqs3test GMP_INCDIR="$GMP_PREFIX/include" \
		GMP_LIB="$GMP_PREFIX/lib/libgmp.a" DETECTED_OS=Linux >/dev/null 2>&1
}

if [ "$KMPQS_ONLY" != 1 ]; then
	build_lasieve lasieve_asm 0
	build_lasieve lasieve_noasm 1
fi

# --- compare ----------------------------------------------------------------
cd "$WORK"
counters() { grep -E '^(success|Stat|  ff)' | sed -E 's/ \([0-9]+\)$//'; }
secs() { grep -E 'total:|^time:' | grep -oE '[0-9]+\.[0-9]+ ?s' | head -1 | tr -d ' s'; }
best_of() {    # best_of cmd... : fastest of $REPS runs
	local r t best=
	for r in $(seq ${REPS:-5}); do
		t=$("$@" | secs)
		if [ -z "$best" ] || awk -v a="$t" -v b="$best" 'BEGIN { exit !(a < b) }'; then best=$t; fi
	done
	echo "$best"
}

status=0
printf "%5s %6s  %-8s %10s %10s %7s\n" bits count counters lasieve_s kmpqs_s ratio
for b in $SIZES; do
	n=${COUNT:-$([ "$b" -gt 100 ] && echo 500 || echo 2000)}
	in="in_${b}.txt"
	./kmpqs_bench -gen "$b" "$n" -dump > "$in"
	if [ "$b" -le 96 ]; then v=-1; t=mpqstest; else v=-3; t=mpqs3test; fi

	k=$(./kmpqs_bench_stat $v "$in")
	if echo "$k" | grep -q '^wrong: [1-9]'; then
		echo "$b: kmpqs returned a wrong factorization"; status=1
	fi
	kt=$(./kmpqs_bench $v -reps ${REPS:-5} "$in" | secs)

	if [ "$KMPQS_ONLY" = 1 ]; then
		cmp="-"; lt="-"; ratio="-"
	else
		if [ "$b" -le 100 ]; then
			if [ "$(./lasieve_noasm/$t "$in" | counters)" = "$(echo "$k" | counters)" ]; then
				cmp=match
			else
				cmp=DIFF; status=1
			fi
			lt=$(best_of ./lasieve_asm/$t "$in")
			ratio=$(awk -v a="$lt" -v k="$kt" 'BEGIN { printf "%.2f", k / a }')
		else
			# lasieve's mpqs3 mostly fails here (and takes seconds per input)
			cmp=skip; lt="-"; ratio="-"
		fi
	fi
	printf "%5s %6s  %-8s %10s %10s %7s  %s\n" "$b" "$n" "$cmp" "$lt" "$kt" "$ratio" \
		"$(echo "$k" | head -1)"
done
exit $status
