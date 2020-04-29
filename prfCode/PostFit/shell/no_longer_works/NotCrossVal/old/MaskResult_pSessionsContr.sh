#!/bin/bash

# $pVOX the minimum proportion of sessions that needed to contribute
# $RTH the R2 cap so bash knows which files to look at
# $3 the subfolder in Results where we'll look for data
# $4 monkey

pVOX=$1
RTH=$2
MONKEY=$3

## old way
# FLD=$3
# MONKEY=$4

# NB The averaged result for this R2 threshold must exist
# Run one of the AverageResults scripts first if it doesn't

# e.g., MaskResult_pSessionsContr 0.2 5 LISA danny 

startpath=$(pwd)

cd ../Results/$MONKEY/AveragedResults_SD/Thr_$RTH
# cd ../Results/$FLD/$MONKEY/AveragedResults/Thr_$RTH

echo 'Only include voxels for which '$pVOX 'sessions contributed (R2_min='$RTH')'

fslmaths pVox_th$RTH.nii.gz -thr $pVOX -bin pvox_mask_$pVOX.nii.gz 
fslmaths pvox_mask_$pVOX.nii.gz -sub 1 -mul 99 fix_mask_$pVOX.nii.gz 

fslmaths MeanAngle_th$RTH.nii.gz -mas pvox_mask_$pVOX.nii.gz -add fix_mask_$pVOX.nii.gz MeanAngle_th${RTH}_pv$pVOX.nii.gz 
fslmaths MeanEccentricity_th$RTH.nii.gz -mas pvox_mask_$pVOX.nii.gz -add fix_mask_$pVOX.nii.gz MeanEccentricity_th${RTH}_pv$pVOX.nii.gz 
fslmaths MeanRFS_th$RTH.nii.gz -mas pvox_mask_$pVOX.nii.gz -add fix_mask_$pVOX.nii.gz MeanRFS_th${RTH}_pv$pVOX.nii.gz 

fslmaths StdAngle_th$RTH.nii.gz -mas pvox_mask_$pVOX.nii.gz -add fix_mask_$pVOX.nii.gz StdAngle_th${RTH}_pv$pVOX.nii.gz 
fslmaths StdEccentricity_th$RTH.nii.gz -mas pvox_mask_$pVOX.nii.gz -add fix_mask_$pVOX.nii.gz StdEccentricity_th${RTH}_pv$pVOX.nii.gz 
fslmaths StdRFS_th$RTH.nii.gz -mas pvox_mask_$pVOX.nii.gz -add fix_mask_$pVOX.nii.gz StdRFS_th${RTH}_pv$pVOX.nii.gz 

cd $startpath