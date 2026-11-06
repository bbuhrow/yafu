#!/bin/bash
# Time one SIQS job on an aarch64 machine under several kernel configs:
# portable (-noneon), and NEON with each med sieve variant (YAFU_NEON_MED).
# Same seed, 1 thread, -noopt.  Uses perf stat cycles:u when perf works,
# else wall time only.  Run from the yafu directory.
#
# Usage: tools/neon_bench.sh '<expr>' [reps] [extra yafu args]
#   e.g. tools/neon_bench.sh 'siqs(rsa(240))' 2
#        tools/neon_bench.sh 'siqs(rsa(240))' 2 -siqsMid sparse
#        CONFIGS="portable split" tools/neon_bench.sh 'siqs(rsa(240))'
expr="$1"; reps=${2:-1}
shift; [ $# -gt 0 ] && shift
extra="$*"
# noneon: every kernel portable; portable/ranges/split: NEON with that med
# sieve (YAFU_NEON_MED); bfixed/bloop: NEON, mid-range buckets
# (-siqsMid bucket) with that NEON bucket fill (YAFU_NEON_BUCKETS)
configs=${CONFIGS:-"noneon portable ranges split"}
perf_ok=0
perf stat -e cycles:u true > /dev/null 2>&1 && perf_ok=1
# yafu ignores the command-line expression when stdin is not a tty, so
# always feed it on stdin
exprf=$(mktemp)
echo "$expr" > "$exprf"
trap 'rm -f "$exprf"' EXIT

printf "%-8s %4s %10s %12s %10s\n" config rep "QS s" "Mcyc" polys
for cfg in $configs; do
    for r in $(seq 1 "$reps"); do
        case $cfg in
            noneon) flag="-noneon"; envv="" ;;
            # mid-range buckets with the NEON bucket fill, either fill mode
            bfixed) flag="-siqsMid bucket"; envv="YAFU_NEON_BUCKETS=fixed" ;;
            bloop)  flag="-siqsMid bucket"; envv="YAFU_NEON_BUCKETS=loop" ;;
            *)      flag=""; envv="YAFU_NEON_MED=$cfg" ;;
        esac
        rm -f siqs.dat
        out=$(mktemp)
        if [ $perf_ok = 1 ]; then
            env $envv perf stat -x, -e cycles:u -o "$out.perf" \
                ./yafu -threads 1 -seed 1 -noopt -v $flag $extra < "$exprf" > "$out" 2>&1
            cyc=$(awk -F, '/cycles/ { printf "%.0f", $1 / 1e6 }' "$out.perf")
        else
            env $envv ./yafu -threads 1 -seed 1 -noopt -v $flag $extra < "$exprf" > "$out" 2>&1
            cyc="-"
        fi
        qs=$(grep -m1 "QS elapsed time" "$out" | sed 's/.*= //; s/ sec.*//')
        polys=$(grep -m1 "total polynomials" "$out" | sed 's/[^0-9]*\([0-9]*\) total.*/\1/')
        printf "%-8s %4s %10s %12s %10s\n" $cfg $r "$qs" "$cyc" "$polys"
        rm -f "$out" "$out.perf"
    done
done
