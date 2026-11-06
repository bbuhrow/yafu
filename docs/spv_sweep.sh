#!/bin/bash
# =============================================================================
# docs/spv_sweep.sh — SIQS small-prime-variation (SPV) tuning sweep
#
# Re-measurement harness for the experiment documented in
# docs/siqs_spv_plan.md (Commits 1-4).  Measures, per (input, small_limit,
# TFSm offset), the total SIQS time, stage timers and relation count, using
# fresh savefiles for every run (a savefile restart silently disables the
# online tuner, and -noopt is used anyway so the formula tf_small_cutoff
# applies).
#
# Usage:
#   cd <scratch dir containing nothing precious>   # runs create ./siqs.dat
#   docs/spv_sweep.sh <path-to-siqs_demo> <input-file> \
#       [limits] [offsets] [extra args...]
#
#   <input-file>      text file containing the number to factor
#   [limits]          SPV limits to sweep           (default "256 512 1024")
#   [offsets]         TFSm offsets vs formula value (default "-10 -5 0 5 10")
#   [extra args...]   passed through to siqs_demo, e.g. -threads 4 -v
#
# Per limit the script first runs a probe (-noopt, no TFSm) and parses the
# guarded "SPV formula tf_small_cutoff = N" diagnostic (requires -v -v, which
# the script adds), then sweeps -siqsTFSm <formula+offset> around it.
#
# Output: one CSV line per run on stdout and in spv_sweep_results.csv:
#   input,limit,tfsm,total,tdiv,msieve,polyupdate,rels,factors_ok
#
# Notes:
# - Build with stage timers if you want the tdiv/msieve/polyupdate columns
#   filled:  make siqs USER_CFLAGS=-DSIQS_PROFILE   (SIQS.c accumulators and
#   the timer block are gated on SIQS_PROFILE; without it those columns are
#   empty and only total/rels are meaningful).
# - The DLP co-tuning variant measured in the plan (subtracting
#   tf_small_cutoff from closnuf in the DLP branch of siqs_static_init) is a
#   commented-out guarded line in factor/qs/SIQS.c; uncomment, rebuild, and
#   re-run this script to re-measure it on other hardware.
# - Single-thread runs give the most stable comparisons; multi-thread runs
#   work but add scheduler noise.
# =============================================================================

SIQS=$1
INFILE=$2
LIMITS=${3:-"256 512 1024"}
OFFSETS=${4:-"-10 -5 0 5 10"}
if [ $# -gt 4 ]; then
    shift 4
    EXTRA="$@"
else
    EXTRA=""
fi

if [ ! -x "$SIQS" ] || [ ! -f "$INFILE" ]; then
    echo "usage: $0 <siqs_demo binary> <input file> [limits] [offsets] [extra args]" >&2
    exit 1
fi

N=$(cat "$INFILE")
CSV=spv_sweep_results.csv
FCUTFILE=$(mktemp)
trap 'rm -f "$FCUTFILE"' EXIT
echo "input,limit,tfsm,total,tdiv,msieve,polyupdate,rels,factors_ok" | tee -a "$CSV"

run_one () {
    # $1 = limit, $2 = tfsm ("" for probe)
    local LIMIT=$1 TFSM=$2 TAG
    rm -f siqs.dat
    # options must be separate argv entries; the expression is separate
    local -a CMD=("$SIQS" "siqs(@)" -siqsSPV "$LIMIT" -noopt -v -v $EXTRA)
    [ -n "$TFSM" ] && CMD+=(-siqsTFSm "$TFSM")
    local OUT
    OUT=$(echo "$N" | "${CMD[@]}" 2>&1)
    local CUT TOTAL TDIV MS POLY RELS OK
    CUT=$(echo "$OUT" | grep -o "trial factoring cutoff at [0-9]*" | grep -o "[0-9]*$")
    echo "$OUT" | grep -o "formula tf_small_cutoff = [0-9]*" | grep -o "[0-9]*$" > "$FCUTFILE"
    TOTAL=$(echo "$OUT" | grep "SIQS elapsed" | grep -oE "[0-9]+\.[0-9]+" | head -1)
    TDIV=$(echo "$OUT" | grep "Tdiv time" | grep -oE "[0-9]+\.[0-9]+")
    MS=$(echo "$OUT" | grep "Med sieve" | grep -oE "[0-9]+\.[0-9]+")
    POLY=$(echo "$OUT" | grep "Poly update" | grep -oE "[0-9]+\.[0-9]+")
    RELS=$(echo "$OUT" | grep -o "[0-9]* rels found" | head -1 | grep -o "^[0-9]*")
    OK=$(echo "$OUT" | grep -c "factors found")
    TAG=$(basename "$INFILE" .txt)
    echo "$TAG,$LIMIT,${TFSM:--},$TOTAL,$TDIV,$MS,$POLY,$RELS,$OK" | tee -a "$CSV"
}

for L in $LIMITS; do
    # probe: discover the formula tf_small_cutoff for this (input, limit)
    PROBE=$(run_one "$L" "")
    FCUT=$(cat "$FCUTFILE")
    echo "# probe $PROBE" >&2
    if [ -z "$FCUT" ]; then
        echo "warning: no 'SPV formula tf_small_cutoff' in probe output for limit $L;" >&2
        echo "         is TRY_SMALLP_OPT defined and the binary built with the diagnostic?" >&2
        continue
    fi
    for OFF in $OFFSETS; do
        T=$((FCUT + OFF))
        run_one "$L" "$T"
    done
done
