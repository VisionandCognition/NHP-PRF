#!/bin/bash
MONKEY=$1
datafld=../Results/$MONKEY/Reg2Anat/BestVoxEstimate

echo Thresholding best-fit pRF-maps for ${MONKEY}

for TH in {1..10}; do
    echo "Threshold R2 > ${TH}" 
    mkdir -p ${datafld}/Thr_${TH}
    fslmaths ${datafld}/R2.nii.gz -thr ${TH} -bin ${datafld}/Thr_${TH}/R2_mask.nii.gz
    fslmaths ${datafld}/Thr_${TH}/R2_mask.nii.gz -sub 1 -mul -1  ${datafld}/Thr_${TH}/R2_negmask.nii.gz
    
    fslmaths ${datafld}/ANG.nii.gz -mas ${datafld}/Thr_${TH}/R2_mask.nii.gz \
        ${datafld}/Thr_${TH}/ANG_th${TH}.nii.gz
    fslmaths ${datafld}/Thr_${TH}/R2_negmask.nii.gz -mul -99 -add ${datafld}/Thr_${TH}/ANG_th${TH}.nii.gz  \
        ${datafld}/Thr_${TH}/ANG_th${TH}.nii.gz
    
    fslmaths ${datafld}/ECC.nii.gz -mas ${datafld}/Thr_${TH}/R2_mask.nii.gz \
        ${datafld}/Thr_${TH}/ECC_th${TH}.nii.gz
    fslmaths ${datafld}/Thr_${TH}/R2_negmask.nii.gz -mul -99 -add ${datafld}/Thr_${TH}/ECC_th${TH}.nii.gz  \
        ${datafld}/Thr_${TH}/ECC_th${TH}.nii.gz

    fslmaths ${datafld}/RFS.nii.gz -mas ${datafld}/Thr_${TH}/R2_mask.nii.gz \
        ${datafld}/Thr_${TH}/RFS_th${TH}.nii.gz
    fslmaths ${datafld}/Thr_${TH}/R2_negmask.nii.gz -mul -99 -add ${datafld}/Thr_${TH}/RFS_th${TH}.nii.gz  \
        ${datafld}/Thr_${TH}/RFS_th${TH}.nii.gz

    fslmaths ${datafld}/SD.nii.gz -mas ${datafld}/Thr_${TH}/R2_mask.nii.gz \
        ${datafld}/Thr_${TH}/SD_th${TH}.nii.gz
    fslmaths ${datafld}/Thr_${TH}/R2_negmask.nii.gz -mul -99 -add ${datafld}/Thr_${TH}/SD_th${TH}.nii.gz  \
        ${datafld}/Thr_${TH}/SD_th${TH}.nii.gz
done

