#!/bin/bash

# $1 the R2 cap so bash knows which files to look at
# $2 the subfolder in Results where we'll look for data
# $3 the monkey name

RTH=$1
MONKEY=$2

## this was the old way with results subfolders. Everything has now been moved to Results root
# FLD=$2
# MONKEY=$3

# The script creates R2-masked versions per session first 

# e.g., CatSessResults.sh 5 us_motreg danny 
#       CatSessResults.sh 5 raw danny

startpath=$(pwd)
cd ../Results/$MONKEY
# cd ../Results/$FLD/$MONKEY

echo "Masking results at R2 value of" $RTH

for sess in 20* ; do
    cd $sess
    
    # if on a mac convert 64bit to 32bit nifti
    if [[ $(uname -s) == Darwin ]]; then
    	mri_convert 'Sess-'$sess'_R2.nii.gz' 'Sess-'$sess'_R2.nii.gz'
    fi
    # create a mask volume based on R2
    fslmaths 'Sess-'$sess'_R2.nii.gz' -nan -thr $RTH -bin 'mask_th'$RTH'.nii.gz'

    # mask results
    fslmaths 'Sess-'$sess'_R2.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' 'R2_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_ang.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' 'ANG_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_ecc.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 'ECC_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_rfsize.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 'RFS_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_sd.nii' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 'SD_th'$RTH'.nii.gz'
    cd ..
done

mkdir CatSessions
cd CatSessions
echo "Concatenating sessions"

res1='SessANG_th'$RTH'.nii.gz'
fslmerge -t $res1 ../*/ANG_th${RTH}.nii.gz

res2='SessECC_th'$RTH'.nii.gz'
fslmerge -t $res2 ../*/ECC_th${RTH}.nii.gz

res3='SessRFS_th'$RTH'.nii.gz'
fslmerge -t $res3 ../*/RFS_th${RTH}.nii.gz

res4='SessR2_th'$RTH'.nii.gz'
fslmerge -t $res4 ../*/R2_th${RTH}.nii.gz

res4='SessSD_th'$RTH'.nii.gz'
fslmerge -t $res4 ../*/SD_th${RTH}.nii.gz

#fsleyes ./REF/HiRes_inResultSpace.nii.gz 'SessANG_th'$RTH'.nii.gz' &
cd $startpath