#!/bin/bash
# Compare yafu's stage-1 / stage-2 lanes with the exact oracle (bs_oracle.py).
# Suyama curves (-avxecm_curve 0; the default is now the Edwards curves).
# usage: oracle_cmp.sh yafu seed pbits nbits B1 B2 nbatch S [extra yafu args]
#   S: 1 plain; S > 1 power x^S; S < -1 Dickson (passed as -avxecm_bs S)
Y=$1; SEED=$2; PB=$3; NB=$4; B1=$5; B2=$6; NBATCH=$7; S=$8; shift 8
read P N <<< "$(python3 mkn.py $SEED $PB $NB)"
BSARG=""
if [ "$S" -lt -1 ]; then BSARG="-avxecm_bs d$((-S))"; elif [ "$S" -gt 1 ]; then BSARG="-avxecm_bs $S"; fi
mism=0; tot1=0; tot2=0
for ((b = 0; b < NBATCH; b++)); do
    s=$((1000 + 8 * b))
    out=$(echo "ecm($N,8)" | $Y -B1ecm $B1 -B2ecm $B2 -sigma $s -v -v -threads 1 -avxecm_stg2 1 -avxecm_curve 0 $BSARG "$@" 2>&1)
    cfg=$(echo "$out" | grep -o "FFT stage 2: B2 = [0-9]*, D = [0-9]*, dF = [0-9]*.*block" | head -1)
    D=$(echo "$cfg" | sed 's/.*D = \([0-9]*\).*/\1/'); DF=$(echo "$cfg" | sed 's/.*dF = \([0-9]*\).*/\1/')
    BL=$(echo "$cfg" | sed 's/.*baby roots), \([0-9]*\) block.*/\1/')
    y1=$(echo "$out" | grep -o "in stage 1 in thread [0-9]*, vec position [0-9]*" | awk '{print $NF}' | sort -n | tr '\n' ' ')
    y2=$(echo "$out" | grep -o "in stage 2 in thread [0-9]*, vec position [0-9]*" | awk '{print $NF}' | sort -n | tr '\n' ' ')
    o=$(python3 bs_oracle.py lanes $P $B1 $D $DF $BL $S $s 1)
    o1=$(echo "$o" | sed 's/.*stage1 \(.*\) stage2.*/\1/' | tr -s ' ' | sed 's/^ //; s/ $//')
    o2=$(echo "$o" | sed 's/.*stage2 *//' | tr -s ' ' | sed 's/^ //; s/ $//')
    y1=$(echo "$y1" | sed 's/ $//'); y2=$(echo "$y2" | sed 's/ $//')
    tot1=$((tot1 + $(echo $y1 | wc -w))); tot2=$((tot2 + $(echo $y2 | wc -w)))
    if [ "$y1" != "$o1" ] || [ "$y2" != "$o2" ]; then
        mism=$((mism + 1))
        echo "  MISMATCH sigma0 $s (D $D dF $DF blocks $BL): yafu s1[$y1] s2[$y2]  oracle s1[$o1] s2[$o2]"
    fi
done
echo "S=$S: $NBATCH batches, yafu stage-1 finds $tot1, stage-2 finds $tot2, batches differing from the oracle: $mism"
