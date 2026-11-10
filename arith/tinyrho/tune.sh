#!/bin/bash
# U / G sweep of the x8 lists (IFMA build): tune.sh [num]
set -e
cd /mnt/c/Users/bbuhr/projects/factoring/yafu-3.1.8/arith/tinyrho
R=$(realpath ../..)
G=$(realpath $R/../gmp-install/6.2.1)
D="-O2 -march=native -DUSE_AVX2 -DUSE_BMI2 -DUSE_AVX512F -DIFMA -I$R/factor/gmp-ecm -I$G/include"
num=${1:-1000}
mkdir -p tune
for cfg in ${CFGS:-"OFF=3" "OFF=2" "OFF=1" "OFF=0" "OFF=2,MAX=12" "OFF=1,MAX=12" "OFF=1,MIN=7" "U104=1" "U104=2" "U104=3" "U104=4"}; do
    F=""; k=""
    for kv in ${cfg//,/ }; do
        k=${kv%=*}; v=${kv#*=}
        case $k in
            U52) F="$F -DTRHO52_U=$v";; U104) F="$F -DTRHO104_U=$v";;
            MIN) F="$F -DTRHO_VB_MIN=$v";; OFF) F="$F -DTRHO_VB_OFF=$v";; MAX) F="$F -DTRHO_VB_MAX=$v";;
            SOFF) F="$F -DTRHO_SB_OFF=$v";; SMAX) F="$F -DTRHO_SB_MAX=$v";;
        esac
    done
    gcc $D $F -o tune/rt rho_test.c $R/factor/gmp-ecm/tinyrho.c $G/lib/libgmp.a -lm
    echo "== $cfg"
    cd $R
    case $k in
        U104) arith/tinyrho/tune/rt u $((num / 2)) 16 28 100 | awk 'NR>1{printf "%s:%s ", $3, $8} END{print ""}';;
        SOFF|SMAX) arith/tinyrho/tune/rt f64 $((num / 2)) | awk 'NR>1{printf "%s:%s ", $1, $5} END{print ""}';;
        *)    arith/tinyrho/tune/rt f52 $num | awk 'NR>1{printf "%s:%s ", $1, $8} END{print ""}';;
    esac
    cd - >/dev/null
done
