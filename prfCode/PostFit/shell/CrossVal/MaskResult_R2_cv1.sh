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
cd ../Results/$MONKEY/$sess
brainmask=../../Reference/${MONKEY}/input_files/Func/brainmask.nii.gz

echo "Masking results for" $MONKEY $sess "at mean xval R2 value of" $RTH

# Create the appropriate masks
fslmaths ${brainmask} -bin brainmask.nii.gz
fslmaths brainmask.nii.gz -sub 1 -mul -1 nonbrainmask.nii.gz
fslmaths nonbrainmask.nii.gz -mul -99 nanreplacement.nii.gz
   
# create average/max/min R2 values
fslmaths 'Sess-'$sess'_R2_1.nii.gz' -nan 'Sess-'$sess'_R2_1.nii.gz'
fslmaths 'Sess-'$sess'_R2_2.nii.gz' -nan 'Sess-'$sess'_R2_2.nii.gz'
fslmaths 'Sess-'$sess'_R2_1.nii.gz' -min 'Sess-'$sess'_R2_2.nii.gz' 'Sess-'$sess'_minR2.nii.gz'
fslmaths 'Sess-'$sess'_R2_1.nii.gz' -max 'Sess-'$sess'_R2_2.nii.gz' 'Sess-'$sess'_maxR2.nii.gz'
fslmaths 'Sess-'$sess'_R2_1.nii.gz' -add 'Sess-'$sess'_R2_2.nii.gz' -div 2 'Sess-'$sess'_meanR2.nii.gz'

# create a mask volume based on mean R2
fslmaths 'Sess-'$sess'_meanR2.nii.gz' -thr $RTH -bin 'mask_th'$RTH'.nii.gz'
fslmaths 'mask_th'$RTH'.nii.gz' -sub 1 -mul -1 'invmask_th'$RTH'.nii.gz'

# mask results
for part in 1 2; do
    echo 'Part ' ${part}
    fslmaths 'Sess-'$sess'_ang_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'ANG_'${part}'_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_ecc_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'ECC_'${part}'_th'$RTH'.nii.gz'
    if [ ${sess:0:3} = css ]; then
        fslmaths 'Sess-'$sess'_expt_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'EXPT_'${part}'_th'$RTH'.nii.gz'
        fslmaths 'Sess-'$sess'_exptsq_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'EXPTSQ_'${part}'_th'$RTH'.nii.gz'
    fi
    fslmaths 'Sess-'$sess'_imag_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'IMAG_'${part}'_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_real_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'REAL_'${part}'_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_FWHM_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'FWHM_'${part}'_th'$RTH'.nii.gz'
    fslmaths 'Sess-'$sess'_rfsize_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'RFS_'${part}'_th'$RTH'.nii.gz'
    if [ ${sess:0:3} = dog ]; then
        fslmaths 'Sess-'$sess'_normamp_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'NAMP_'${part}'_th'$RTH'.nii.gz'
        fslmaths 'Sess-'$sess'_sd2_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -div 10 -add nanreplacement.nii.gz 'IRFS_'${part}'_th'$RTH'.nii.gz'
        fslmaths 'Sess-'$sess'_sdratio_'${part}'.nii.gz' -nan -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'SDRATIO_'${part}'_th'$RTH'.nii.gz'
    fi

    fslmaths 'Sess-'$sess'_ang_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_ang_'${part}'.nii.gz'
    fslmaths 'Sess-'$sess'_ecc_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_ecc_'${part}'.nii.gz'
    if [ ${sess:0:3} = css ]; then
        fslmaths 'Sess-'$sess'_expt_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_expt_'${part}'.nii.gz' 
        fslmaths 'Sess-'$sess'_exptsq_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_exptsq_'${part}'.nii.gz' 
    fi
    fslmaths 'Sess-'$sess'_imag_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_imag_'${part}'.nii.gz' 
    fslmaths 'Sess-'$sess'_real_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_real_'${part}'.nii.gz' 
    fslmaths 'Sess-'$sess'_FWHM_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_FWHM_'${part}'.nii.gz'  
    fslmaths 'Sess-'$sess'_rfsize_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_rfsize_'${part}'.nii.gz'
    if [ ${sess:0:3} = dog ]; then
        fslmaths 'Sess-'$sess'_normamp_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_normamp_'${part}'.nii.gz'
        fslmaths 'Sess-'$sess'_sd2_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_sd2_'${part}'.nii.gz'
        fslmaths 'Sess-'$sess'_sdratio_'${part}'.nii.gz' -nan -add nanreplacement.nii.gz 'Sess-'$sess'_sdratio_'${part}'.nii.gz'
    fi

    mkdir -p TH_$RTH
    mv 'ANG_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    mv 'ECC_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    if [ ${sess:0:3} = css ]; then
        mv 'EXPT_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
        mv 'EXPTSQ_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    fi
    mv 'IMAG_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    mv 'REAL_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    mv 'FWHM_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    mv 'RFS_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    if [ ${sess:0:3} = dog ]; then
        mv 'NAMP_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
        mv 'IRFS_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
        mv 'SDRATIO_'${part}'_th'$RTH'.nii.gz' ./TH_$RTH/
    fi    
done
mv 'mask_th'$RTH'.nii.gz' ./TH_$RTH/
mv 'invmask_th'$RTH'.nii.gz' ./TH_$RTH/
mv nanreplacement.nii.gz ./TH_$RTH/

cd ./TH_$RTH
# average the two xval results
#fslmaths 'ECC_1_th'$RTH'.nii.gz' -add 'ECC_2_th'$RTH'.nii.gz' -div 2 'ECC_th'$RTH'.nii.gz'
#fslmaths 'IMAG_1_th'$RTH'.nii.gz' -add 'IMAG_2_th'$RTH'.nii.gz' -div 2 'IMAG_th'$RTH'.nii.gz'
#fslmaths 'REAL_1_th'$RTH'.nii.gz' -add 'REAL_2_th'$RTH'.nii.gz' -div 2 'REAL_th'$RTH'.nii.gz'
fslmaths 'FWHM_1_th'$RTH'.nii.gz' -add 'FWHM_2_th'$RTH'.nii.gz' -div 2 'FWHM_th'$RTH'.nii.gz'
fslmaths 'RFS_1_th'$RTH'.nii.gz' -add 'RFS_2_th'$RTH'.nii.gz' -div 2 'RFS_th'$RTH'.nii.gz'
if [ ${sess:0:3} = css ]; then
    fslmaths 'EXPT_1_th'$RTH'.nii.gz' -add 'EXPT_2_th'$RTH'.nii.gz' -div 2 'EXPT_th'$RTH'.nii.gz'
    fslmaths 'EXPTSQ_1_th'$RTH'.nii.gz' -add 'EXPTSQ_2_th'$RTH'.nii.gz' -div 2 'EXPTSQ_th'$RTH'.nii.gz'
fi
if [ ${sess:0:3} = dog ]; then
    fslmaths 'NAMP_1_th'$RTH'.nii.gz' -add 'NAMP_2_th'$RTH'.nii.gz' -div 2 'NAMP_th'$RTH'.nii.gz'
    fslmaths 'IRFS_1_th'$RTH'.nii.gz' -add 'IRFS_2_th'$RTH'.nii.gz' -div 2 'IRFS_th'$RTH'.nii.gz'
    fslmaths 'SDRATIO_1_th'$RTH'.nii.gz' -add 'SDRATIO_2_th'$RTH'.nii.gz' -div 2 'SDRATIO_th'$RTH'.nii.gz'
fi

for xv in 1 2; do
    #fslmaths 'ANG_'${xv}'_th'$RTH'.nii.gz' -mul 0.0174533 -cos -mul -1 -mul 'ECC_'${xv}'_th'$RTH'.nii.gz' 'X_'${xv}'_th'$RTH'.nii.gz'
    #fslmaths 'ANG_'${xv}'_th'$RTH'.nii.gz' -mul 0.0174533 -sin -mul 'ECC_'${xv}'_th'$RTH'.nii.gz' 'Y_'${xv}'_th'$RTH'.nii.gz'
    fslmaths 'ANG_'${xv}'_th'$RTH'.nii.gz' -mul 0.0174533 -cos -mul 'ECC_'${xv}'_th'$RTH'.nii.gz' 'X_'${xv}'_th'$RTH'.nii.gz'
    fslmaths 'ANG_'${xv}'_th'$RTH'.nii.gz' -mul 0.0174533 -sin -mul 'ECC_'${xv}'_th'$RTH'.nii.gz' 'Y_'${xv}'_th'$RTH'.nii.gz'
done

fslmaths 'X_1_th'$RTH'.nii.gz' -add 'X_2_th'$RTH'.nii.gz' -div 2 -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'X_th'$RTH'.nii.gz'
fslmaths 'Y_1_th'$RTH'.nii.gz' -add 'Y_2_th'$RTH'.nii.gz' -div 2 -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'Y_th'$RTH'.nii.gz'

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

fslmaths 'Y_th'$RTH'.nii.gz' -div 'X_th'$RTH'.nii.gz' -atan -mul 57.2957795 'rawANG.nii.gz'
# +X +Y  >>  atan(Y/X)
# do nothing
# -X -Y  >>  atan(Y/X) -180
fslmaths 'negXnegY.nii.gz' -mul -180 -add 'rawANG.nii.gz' 'raw2ANG.nii.gz'
# -X +Y  >>  atan(Y/X) +180
fslmaths 'negXposY.nii.gz' -mul 180 -add 'raw2ANG.nii.gz' -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'raw3ANG.nii.gz'
fslmaths 'raw3ANG.nii.gz' -binv -mul 360 -add 'raw3ANG.nii.gz' -mas 'mask_th'$RTH'.nii.gz' -add nanreplacement.nii.gz 'ANG_th'$RTH'.nii.gz'

# clean up, move some intermediate results to a separate folder
mkdir -p 'helper_vols'
mv *_1_*.nii.gz ./helper_vols/
mv *_2_*.nii.gz ./helper_vols/
mv *pos*.nii.gz ./helper_vols/
mv *neg*.nii.gz ./helper_vols/
mv *raw*.nii.gz ./helper_vols/
mv *X2*.nii.gz ./helper_vols/
mv *Y2*.nii.gz ./helper_vols/

cd $startpath