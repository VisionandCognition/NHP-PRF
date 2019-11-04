#!/bin/bash
MONKEY=eddy

declare -a sessions=(
    'csshrf_cv0_dhrf'
    'csshrf_cv0_mhrf'
    'doghrf_cv0_dhrf'
    'doghrf_cv0_mhrf'
    'linhrf_cv0_dhrf'
    'linhrf_cv0_mhrf'
    'linhrf_cv0_dhrf_neggain'
    'linhrf_cv0_mhrf_neggain'
)

for sess in "${sessions[@]}"; do
    for rth in 0 1 2 4 5 10; do
        ./MaskResult_R2_cv0.sh $MONKEY $sess $rth
    done
done
