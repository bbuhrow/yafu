#!/bin/bash
# SIQS benchmark / correctness matrix.
#
# usage: tools/bench_matrix.sh [options] variant [variant ...]
#   variant: label=path/to/yafu[,extra yafu args]
#     e.g.  avx512=~/yafu_matrix/avx512/yafu
#           a512_sink=~/yafu_matrix/avx512/yafu,-siqsKernel sink
# options:
#   -i "c77 c79"   inputs, by name from tools/bench_inputs.txt (default c77 c79)
#   -n "2 4 8"     NB values; "def" = yafu's choice    (default "4")
#   -r 2           repetitions                           (default 1)
#   -t 1           threads                               (default 1)
#   -s 777         rng seed; "none" = no -seed/-noopt    (default 777)
#   -o file.csv    output CSV                            (default bench_<date>.csv)
#   -w dir         work dir (siqs.dat etc. go here)      (default ~/siqsbench)
#   -p             wrap runs in perf stat (cycles, instructions, branch-misses)
#
# Runs one job at a time.  The variant loop is innermost, so drift spreads
# across variants.  siqs.dat is deleted before every run.  With a seed
# (default), -noopt is also passed, so builds that sieve the same prime
# ranges the same way do bit-identical work; summarize with
# tools/bench_summary.py, which also checks the work counters.

set -e
HERE=$(cd "$(dirname "$0")" && pwd)
INPUTS="c77 c79"; NBS="4"; REPS=1; THREADS=1; SEED=777; PERF=0
OUT="bench_$(date +%Y%m%d_%H%M%S).csv"; WORK="$HOME/siqsbench"

while getopts "i:n:r:t:s:o:w:p" opt; do
    case $opt in
        i) INPUTS=$OPTARG ;; n) NBS=$OPTARG ;; r) REPS=$OPTARG ;;
        t) THREADS=$OPTARG ;; s) SEED=$OPTARG ;; o) OUT=$OPTARG ;;
        w) WORK=$OPTARG ;; p) PERF=1 ;;
        *) sed -n 2,20p "$0"; exit 1 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 0 ] && { sed -n 2,20p "$0"; exit 1; }

if pgrep -x yafu > /dev/null; then
    echo "another yafu is running; refusing to benchmark" >&2
    exit 1
fi

OUT=$(realpath -m "$OUT")
mkdir -p "$WORK"
cd "$WORK"
[ -f "$HERE/../yafu.ini" ] && [ ! -f yafu.ini ] && cp "$HERE/../yafu.ini" .

echo "host,cpu,variant,input,nb,rep,qs_sec,polys,touched,blocks,full,partial,factored,cycles,instructions,branch_misses,args" > "$OUT"
CPU=$(grep -m1 "model name" /proc/cpuinfo 2>/dev/null | sed 's/.*: //; s/,/ /g')

for rep in $(seq 1 "$REPS"); do
for name in $INPUTS; do
    line=$(grep "^$name " "$HERE/bench_inputs.txt") || { echo "no input $name" >&2; exit 1; }
    N=$(echo "$line" | awk '{print $2}'); F1=$(echo "$line" | awk '{print $3}')
    for nb in $NBS; do
        for v in "$@"; do
            label=${v%%=*}; rest=${v#*=}
            bin=${rest%%,*}; extra=""
            [ "$rest" != "$bin" ] && extra=${rest#*,}
            bin=${bin/#\~/$HOME}
            args="-threads $THREADS -v"
            [ "$nb" != "def" ] && args="$args -siqsNB $nb"
            [ "$SEED" != "none" ] && args="$args -seed $SEED -noopt"
            args="$args $extra"
            rm -f siqs.dat
            if [ "$PERF" = 1 ]; then
                echo "$N" | perf stat -x, -o perf.txt -e cycles,instructions,branch-misses \
                    "$bin" "siqs(@)" $args > run.txt 2>&1 || true
            else
                echo "$N" | "$bin" "siqs(@)" $args > run.txt 2>&1 || true
                : > perf.txt
            fi
            tr '\r' '\n' < run.txt > run_n.txt
            qs=$(grep -m1 "^QS elapsed time" run_n.txt | awk '{print $5}')
            polys=$(grep -m1 "total polynomials" run_n.txt | awk '{print $3}')
            touched=$(grep -m1 "trial division touched" run_n.txt | awk '{print $4}')
            blocks=$(grep -m1 "total blocks sieved" run_n.txt | awk '{print $5}' | tr -d ,)
            last=$(grep "rels found" run_n.txt | tail -1)
            full=$(echo "$last" | sed -n 's/.*: \([0-9]*\) full.*/\1/p')
            partial=$(echo "$last" | sed -n 's/.* from \([0-9]*\) partial.*/\1/p')
            fac=0; grep -q "$F1" run_n.txt && fac=1
            cyc=$(grep ",cycles" perf.txt | cut -d, -f1)
            ins=$(grep ",instructions" perf.txt | cut -d, -f1)
            brm=$(grep ",branch-misses" perf.txt | cut -d, -f1)
            echo "$(hostname),$CPU,$label,$name,$nb,$rep,$qs,$polys,$touched,$blocks,$full,$partial,$fac,$cyc,$ins,$brm,$extra" >> "$OUT"
            printf "%-14s %-4s NB%-4s rep%s  %8s s  polys %-7s touched %-9s %s\n" \
                "$label" "$name" "$nb" "$rep" "$qs" "$polys" "$touched" \
                "$([ $fac = 1 ] && echo ok || echo NOT-FACTORED)"
            cp run.txt "last_${label}_${name}_${nb}.txt"
        done
    done
done
done
rm -f siqs.dat
echo "results: $OUT"
