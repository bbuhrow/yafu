#!/bin/bash
# Stage-1 timing sweep of the Edwards stage 1 (-avxecm_curve 2) over the
# window width, chunk length and prefetch (AVXECM_EDW_CFG="w,chunk,pf"),
# with the PRAC stage 1 (Suyama, curve 0) as the reference.
# usage (from this directory): edw_sweep.sh yafu label [core]
# Prints: label bits B1 config seconds (min of 2, 8 curves, 1 thread).
Y=$1; LABEL=$2; CORE=${3:-3}
run() {   # N B1 curve cfg
    local best=
    for r in 1 2; do
        t=$(echo "ecm($1,8)" | AVXECM_EDW_CFG="$4" taskset -c $CORE $Y -B1ecm $2 -B2ecm $(($2 + 1)) \
            -sigma 1000 -v -v -threads 1 -avxecm_curve $3 2>&1 | grep 'Stage 1 took' | awk '{print $4}')
        if [ -z "$best" ] || [ $(echo "$t < $best" | bc) = 1 ]; then best=$t; fi
    done
    echo "$best"
}
for spec in "413 1000000" "1037 1000000" "2048 250000"; do
    set -- $spec; NB=$1; B1=$2
    read P N <<< "$(python3 ../vecpoly/ecm_tests/mkn.py 5 200 $NB)"
    echo "$LABEL $NB $B1 prac $(run $N $B1 0 '')"
    for cfg in 5,16384,1 6,16384,1 7,16384,1 8,16384,1 9,16384,1 10,16384,1 11,16384,1 12,16384,1 \
               8,4096,1 8,65536,1 10,4096,1 10,65536,1 8,16384,0 10,16384,0 12,16384,0; do
        echo "$LABEL $NB $B1 edw:$cfg $(run $N $B1 2 $cfg)"
    done
done
