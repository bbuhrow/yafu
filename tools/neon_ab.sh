#!/bin/bash
# A/B one SIQS job on an aarch64 yafu: NEON kernels vs -noneon (portable),
# same seed, 1 thread, and report whether the deterministic statistics
# (polys, reports, survivors, relations, factors) match.  Runs under
# qemu-aarch64-static on an x86 host, or natively on an ARM machine.
#
# Usage: tools/neon_ab.sh <yafu dir> '<expr>' [extra yafu args]
#   e.g. tools/neon_ab.sh ~/neon/yafu 'siqs(rsa(232))'
#
# -noopt is always passed: the SPV tuner picks tf_small_cutoff from
# measured rels/sec, which makes otherwise identical runs differ.
# -siqsMid sparse is passed unless the extra args set -siqsMid: NEON and
# portable kernels default to different mid strategies (sparse / bucket).
cd "$1" || exit 1
expr="$2"; shift 2
case " $* " in
    *" -siqsMid "*) ;;
    *) set -- -siqsMid sparse "$@" ;;
esac
run=""
[ "$(uname -m)" != aarch64 ] && run=qemu-aarch64-static
out=${TMPDIR:-/tmp}
for mode in neon noneon; do
    flag=""; [ $mode = noneon ] && flag="-noneon"
    rm -f siqs.dat
    # yafu reads a non-tty stdin as a batch file, so feed the expression there
    start=$(date +%s.%N)
    echo "$expr" | $run ./yafu -threads 1 -seed 1 -noopt -v $flag "$@" > "$out/ab_$mode.txt" 2>&1
    end=$(date +%s.%N)
    # (not the "rels found" progress lines: those are printed on a timer)
    grep -E "total polynomials|trial division touched|total reports|total blocks|dlp-ecm|^P[0-9]+ = " \
        "$out/ab_$mode.txt" | sed 's/(.*//' > "$out/ab_$mode.stats"
    printf "== %s (%s s): %s\n" $mode "$(awk "BEGIN { printf \"%.1f\", $end - $start }")" \
        "$(grep -E 'strategy:' "$out/ab_$mode.txt" | sed 's/strategy: //')"
    cat "$out/ab_$mode.stats"
done
if cmp -s "$out/ab_neon.stats" "$out/ab_noneon.stats"; then
    echo "MATCH"
else
    echo "DIFFER"; diff "$out/ab_neon.stats" "$out/ab_noneon.stats"
fi
