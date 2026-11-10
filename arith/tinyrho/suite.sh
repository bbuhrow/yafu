#!/bin/bash
# correctness + timing suite for one ISA build: suite.sh [ISA] [num]
D=/mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8/arith/tinyrho
isa=${1:-ifma}; num=${2:-1000}
bash $D/run.sh $isa f52 $num
bash $D/run.sh $isa f64 $num | tail -6
bash $D/run.sh $isa f128 $num
bash $D/run.sh $isa u $((num / 2)) 16 30 100
bash $D/run.sh $isa u $((num / 4)) 16 26 120
