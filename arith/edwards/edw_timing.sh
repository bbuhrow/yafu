#!/bin/bash
# Stage-1 times, PRAC (Suyama, -avxecm_curve 0) vs the Edwards stage 1
# (-avxecm_curve 2, default w / chunk), 8 curves, 1 thread, min of 3.
# usage (from this directory): edw_timing.sh yafu label [core]
Y=$1; LABEL=$2; CORE=${3:-3}
run() {   # N B1 curve
    local best=
    for r in 1 2 3; do
        t=$(echo "ecm($1,8)" | taskset -c $CORE $Y -B1ecm $2 -B2ecm $(($2 + 1)) \
            -sigma 1000 -v -v -threads 1 -avxecm_curve $3 2>&1 | grep 'Stage 1 took' | awk '{print $4}')
        if [ -z "$best" ] || [ $(echo "$t < $best" | bc) = 1 ]; then best=$t; fi
    done
    echo "$best"
}
for spec in "256 1000000" "413 1000000" "700 1000000" "1037 1000000" "1500 250000" "2048 250000" \
            "413 50000" "1037 50000" "413 11000"; do
    set -- $spec; NB=$1; B1=$2
    read P N <<< "$(python3 ../vecpoly/ecm_tests/mkn.py 5 120 $NB)"
    a=$(run $N $B1 0); b=$(run $N $B1 2)
    echo "$LABEL $NB $B1 prac $a edw $b ratio $(echo "scale=3; $a / $b" | bc)"
done
