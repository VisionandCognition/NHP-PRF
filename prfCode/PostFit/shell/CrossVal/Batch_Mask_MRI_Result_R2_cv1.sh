#!/bin/bash
#MONKEY=danny
#MONKEY=eddy

declare -a monkeys=(
    'danny'
    'eddy'
)

declare -a sessions=(
    'csshrf_cv1_dhrf'
    'csshrf_cv1_mhrf'
    'doghrf_cv1_dhrf'
    'doghrf_cv1_mhrf'
    'linhrf_cv1_dhrf'
    'linhrf_cv1_mhrf'
    'linhrf_cv1_dhrf_neggain'
    'linhrf_cv1_mhrf_neggain'
)

for m in "${monkeys[@]}"; do
    for sess in "${sessions[@]}"; do
        for rth in 0 1 2 4 5 10; do
            ./Mask_MRI_Result_R2_cv1.sh $m $sess $rth
        done
    done
done