#!/bin/bash
# usage: bench_s2.sh yafu N B1 [extra args]
# stage 1 / stage 2 times for: pairing (default B2 = 100 B1), FFT at
# B2 = 100 B1, FFT at its default B2 (GMP-ECM's).
Y=$1; N=$2; B1=$3; shift 3
run() {
    echo "ecm($N,8)" | $Y -B1ecm $B1 -sigma 1000 -v -v -threads 1 "$@" 2>&1 | \
        grep -E "Stage 1 took|Stage 2.*took|FFT stage 2:|FFT stage 2 \(thread" | grep -v "Init took" | \
        sed -e 's/^/    /'
}
echo "  pairing, B2 = 100 B1:"; run -avxecm_stg2 0 "$@"
echo "  FFT, B2 = 100 B1:";     run -avxecm_stg2 1 -B2ecm $((100 * B1)) "$@"
echo "  FFT, default B2:";      run -avxecm_stg2 1 "$@"
