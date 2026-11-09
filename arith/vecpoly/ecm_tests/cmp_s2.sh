#!/bin/bash
# Lane-for-lane comparison of AVX-ECM stage 2 methods with fixed sigmas.
# usage: cmp_s2.sh yafu_binary N B1 B2 nbatch [extra yafu args]
# For each batch (8 curves, sigma0 = 1000 + 8*b): runs pairing (-avxecm_stg2 0)
# and FFT (-avxecm_stg2 1) at the same requested B2, prints the stage-2 lanes
# found by each, and flags any pairing find the FFT run missed.
Y=$1; N=$2; B1=$3; B2=$4; NB=$5; shift 5
miss=0; extra=0; tot_p=0; tot_f=0
for ((b = 0; b < NB; b++)); do
    s=$((1000 + 8 * b))
    p=$(echo "ecm($N,8)" | $Y -B1ecm $B1 -B2ecm $B2 -sigma $s -v -v -threads 1 -avxecm_stg2 0 "$@" 2>&1 | \
        grep -o "in stage 2 in thread [0-9]*, vec position [0-9]*" | awk '{print $NF}' | sort -n | tr '\n' ' ')
    f=$(echo "ecm($N,8)" | $Y -B1ecm $B1 -B2ecm $B2 -sigma $s -v -v -threads 1 -avxecm_stg2 1 "$@" 2>&1 | \
        grep -o "in stage 2 in thread [0-9]*, vec position [0-9]*" | awk '{print $NF}' | sort -n | tr '\n' ' ')
    for l in $p; do
        tot_p=$((tot_p + 1))
        if ! echo " $f " | grep -q " $l "; then miss=$((miss + 1)); echo "  MISS: sigma0 $s lane $l (pair: $p / fft: $f)"; fi
    done
    for l in $f; do
        tot_f=$((tot_f + 1))
        if ! echo " $p " | grep -q " $l "; then extra=$((extra + 1)); fi
    done
done
echo "stage-2 finds: pairing $tot_p, FFT $tot_f; pairing finds missed by FFT: $miss; FFT-only finds: $extra"
