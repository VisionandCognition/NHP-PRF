#!/bin/bash
MONKEY=eddy

declare -a sessions=(
    'doghrf_cv1_dhrf'
    'doghrf_cv1_mhrf'
    'linhrf_cv1_dhrf_neggain'
    'linhrf_cv1_mhrf_neggain'
)

declare -a sessions2=(
    'csshrf_cv1_dhrf'
    'csshrf_cv1_mhrf'
    'doghrf_cv1_dhrf'
    'doghrf_cv1_mhrf'
    'linhrf_cv1_dhrf'
    'linhrf_cv1_mhrf'
    'linhrf_cv1_dhrf_neggain'
    'linhrf_cv1_mhrf_neggain'
)

for sess in "${sessions[@]}"; do
    for rth in 0 1 2 4 5 10; do
        ./MaskResult_R2_cv1.sh $MONKEY $sess $rth
    done
done
