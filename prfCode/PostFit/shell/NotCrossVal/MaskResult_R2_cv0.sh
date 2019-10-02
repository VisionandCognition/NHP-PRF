#!/bin/bash

MONKEY=$1
sess=$2
RTH=$3

# debugging ===================
#MONKEY='danny'
#sess='csshrf_cv1_dhrf'
#RTH=2
# /debugging ==================

startpath=$(pwd)
cd ./Results/$MONKEY/$sess
brainmask=../../Reference/${MONKEY}/input_files/Func/brainmask.nii.gz

echo "Masking results for" $MONKEY $sess "at R2 value of" $RTH

# Create the appropriate masks
fslmaths ${brainmask} -bin brainmask.nii.gz
fslmaths brainmask.nii.gz -sub 1 -mul -1 nonbrainmask.nii.gz
fslmaths nonbrainmask.nii.gz -mul -99 nanreplacement.nii.gz
   
# create average/max/min R2 values
fslmaths 'Sess-'$sess'_R2.nii.gz' -nan 'Sess-'$sess'_R2.nii.gz'

# create a mask volume based on R2
fslmaths 'Sess-'$sess'_R2.nii.gz' -thr $RTH -bin 'mask_th'$RTH'.nii.gz'
fslmaths 'mask_th'$RTH'.nii.gz' -sub 1 -mul -1 'invmask_th'$RTH'.nii.gz'

# mask results
fslmaths 'Sess-'$sess'_ang.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'ANG_th'$RTH'.nii.gz'
fslmaths 'Sess-'$sess'_ecc.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'ECC_th'$RTH'.nii.gz'
if [ ${sess:0:3} = css ]; then
    fslmaths 'Sess-'$sess'_expt.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'EXPT_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_exptsq.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'EXPTSQ_th'$RTH'.nii.gz'
fi
fslmaths 'Sess-'$sess'_imag.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'IMAG_th'$RTH'.nii.gz'
fslmaths 'Sess-'$sess'_real.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'REAL_th'$RTH'.nii.gz'
fslmaths 'Sess-'$sess'_FWHM.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'FWHM_th'$RTH'.nii.gz'
fslmaths 'Sess-'$sess'_rfsize.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'RFS_th'$RTH'.nii.gz'
if [ ${sess:0:3} = dog ]; then
    fslmaths 'Sess-'$sess'_normamp.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'NAMP_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_sd2.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'IRFS_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_sdratio.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'SDRATIO_th'$RTH'.nii.gz'
fi

fslmaths 'Sess-'$sess'_ang.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_ang.nii.gz'
fslmaths 'Sess-'$sess'_ecc.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_ecc.nii.gz'
if [ ${sess:0:3} = css ]; then
    fslmaths 'Sess-'$sess'_expt.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_expt.nii.gz' 
    fslmaths 'Sess-'$sess'_exptsq.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_exptsq.nii.gz' 
fi
fslmaths 'Sess-'$sess'_imag.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_imag.nii.gz' 
fslmaths 'Sess-'$sess'_real.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_real.nii.gz' 
fslmaths 'Sess-'$sess'_FWHM.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_FWHM.nii.gz'  
fslmaths 'Sess-'$sess'_rfsize.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_rfsize.nii.gz'
if [ ${sess:0:3} = dog ]; then
    fslmaths 'Sess-'$sess'_normamp.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_normamp.nii.gz'
    fslmaths 'Sess-'$sess'_sd2.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_sd2.nii.gz'
    fslmaths 'Sess-'$sess'_sdratio.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_sdratio.nii.gz'
fi

mkdir -p TH_$RTH
mv 'ANG_th'$RTH'.nii.gz' ./TH_$RTH/
mv 'ECC_th'$RTH'.nii.gz' ./TH_$RTH/
if [ ${sess:0:3} = css ]; then
    mv 'EXPT_th'$RTH'.nii.gz' ./TH_$RTH/
    mv 'EXPTSQ_th'$RTH'.nii.gz' ./TH_$RTH/
fi
mv 'IMAG_th'$RTH'.nii.gz' ./TH_$RTH/
mv 'REAL_th'$RTH'.nii.gz' ./TH_$RTH/
mv 'FWHM_th'$RTH'.nii.gz' ./TH_$RTH/
mv 'RFS_th'$RTH'.nii.gz' ./TH_$RTH/
if [ ${sess:0:3} = dog ]; then
    mv 'NAMP_th'$RTH'.nii.gz' ./TH_$RTH/
    mv 'IRFS_th'$RTH'.nii.gz' ./TH_$RTH/
    mv 'SDRATIO_th'$RTH'.nii.gz' ./TH_$RTH/
fi    

mv 'mask_th'$RTH'.nii.gz' ./TH_$RTH/
mv 'invmask_th'$RTH'.nii.gz' ./TH_$RTH/
mv nanreplacement.nii.gz ./TH_$RTH/

cd ./TH_$RTH

fslmaths 'ANG_th'$RTH'.nii.gz' -mul 0.0174533 -cos -mul 'ECC_th'$RTH'.nii.gz' 'X_th'$RTH'.nii.gz'
fslmaths 'ANG_th'$RTH'.nii.gz' -mul 0.0174533 -sin -mul 'ECC_th'$RTH'.nii.gz' 'Y_th'$RTH'.nii.gz'

fslmaths 'X_th'$RTH'.nii.gz' -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'X_th'$RTH'.nii.gz'
fslmaths 'Y_th'$RTH'.nii.gz' -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'Y_th'$RTH'.nii.gz'

fslmaths 'X_th'$RTH'.nii.gz' -sqr 'X2_th'$RTH'.nii.gz'
fslmaths 'Y_th'$RTH'.nii.gz' -sqr 'Y2_th'$RTH'.nii.gz'
fslmaths 'X2_th'$RTH'.nii.gz' -add 'Y2_th'$RTH'.nii.gz' -sqrt -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'ECC_th'$RTH'.nii.gz'

# calculate angle polar to cartesian
fslmaths 'X_th'$RTH'.nii.gz' -binv -mas 'mask_th'$RTH'.nii.gz' 'negX.nii.gz'
fslmaths 'Y_th'$RTH'.nii.gz' -binv -mas 'mask_th'$RTH'.nii.gz' 'negY.nii.gz'
fslmaths 'X_th'$RTH'.nii.gz' -bin -mas 'mask_th'$RTH'.nii.gz' 'posX.nii.gz'
fslmaths 'Y_th'$RTH'.nii.gz' -bin -mas 'mask_th'$RTH'.nii.gz' 'posY.nii.gz'
fslmaths 'negX.nii.gz' -mul 'posY.nii.gz' -mas 'mask_th'$RTH'.nii.gz' 'negXposY.nii.gz'
fslmaths 'negX.nii.gz' -mul 'negY.nii.gz' -mas 'mask_th'$RTH'.nii.gz' 'negXnegY.nii.gz'

fslmaths 'Y_th'$RTH'.nii.gz' -div 'ECC_th'$RTH'.nii.gz' -atan -mul 57.2957795 'rawANG.nii.gz'
# +X +Y  >>  atan(Y/X)
# do nothing
# -X -Y  >>  atan(Y/X) -180
fslmaths 'negXnegY.nii.gz' -mul -180 -add 'rawANG.nii.gz' 'raw2ANG.nii.gz'
# -X -Y  >>  atan(Y/X) +180
fslmaths 'negXposY.nii.gz' -mul 180 -add 'raw2ANG.nii.gz' -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'raw3ANG.nii.gz'
fslmaths 'raw3ANG.nii.gz' -binv -mul 360 -add 'raw3ANG.nii.gz' -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'ANG_th'$RTH'.nii.gz'

# clean up, move some intermediate results to a separate folder
mkdir -p 'helper_vols'
mv *pos*.nii.gz ./helper_vols/
mv *neg*.nii.gz ./helper_vols/
mv *raw*.nii.gz ./helper_vols/
mv *X2*.nii.gz ./helper_vols/
mv *Y2*.nii.gz ./helper_vols/

cd $startpath