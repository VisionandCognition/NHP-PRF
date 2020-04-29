#!/bin/bash
MONKEY=danny

startpath=$(pwd)
cd ../Results/$MONKEY/
mkdir -p BestVoxEstimate
echo Processing pRF maps from ${MONKEY}

mkdir -p ./BestVoxEstimate/Temp

R2=./CatSessions/SessR2.nii.gz
# make sure this isn't over concatenated
fslroi ${R2} ${R2} 0 16

# get a non-brain mask
#cp /NHP_MRI/NHP-BIDS/manual-masks/final/sub-danny/ses-20180117/func/T1_to_func_brainmask_zcrop.nii.gz \
#    ./BestVoxEstimate/brainmask.nii.gz
fslmaths ./BestVoxEstimate/brainmask.nii.gz -sub 1 -mul -1 ./BestVoxEstimate/non-brainmask.nii.gz

# get the session-concatenated maps 
ANG=./CatSessions/SessANG_th0.nii.gz
ECC=./CatSessions/SessECC_th0.nii.gz
RFS=./CatSessions/SessRFS_th0.nii.gz
SD=./CatSessions/SessSD_th0.nii.gz

# get the max R2 per voxel
fslmaths ${R2} -Tmax ./BestVoxEstimate/R2_max.nii.gz
fslmaths ${R2} -Tmaxn ./BestVoxEstimate/R2_max_idx.nii.gz

# convert the 4D max R2 index volume to session based masks
sess_arr=( 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 )

angstr= 
eccstr= 
rfsstr= 
sdstr=

# split the concatenate results
fslsplit ${ANG} ./BestVoxEstimate/Temp/ANG -t
fslsplit ${ECC} ./BestVoxEstimate/Temp/ECC -t
fslsplit ${RFS} ./BestVoxEstimate/Temp/RFS -t
fslsplit ${SD} ./BestVoxEstimate/Temp/SD -t

for sess in "${sess_arr[@]}"; do
    if [ "${sess}" = "00" ]; then
        fslmaths ./BestVoxEstimate/R2_max_idx.nii.gz -uthr 0.5 -add 1 \
            -bin ./BestVoxEstimate/Temp/R2mask_${sess}.nii.gz
    else
        fslmaths ./BestVoxEstimate/R2_max_idx.nii.gz -thr ${sess} -uthr ${sess} \
            -bin ./BestVoxEstimate/Temp/R2mask_${sess}.nii.gz
    fi

    fslmaths ./BestVoxEstimate/Temp/ANG00${sess}.nii.gz -mas ./BestVoxEstimate/Temp/R2mask_${sess}.nii.gz \
        ./BestVoxEstimate/Temp/mANG00${sess}.nii.gz
    fslmaths ./BestVoxEstimate/Temp/ECC00${sess}.nii.gz -mas ./BestVoxEstimate/Temp/R2mask_${sess}.nii.gz \
        ./BestVoxEstimate/Temp/mECC00${sess}.nii.gz
    fslmaths ./BestVoxEstimate/Temp/RFS00${sess}.nii.gz -mas ./BestVoxEstimate/Temp/R2mask_${sess}.nii.gz \
        ./BestVoxEstimate/Temp/mRFS00${sess}.nii.gz
    fslmaths ./BestVoxEstimate/Temp/SD00${sess}.nii.gz -mas ./BestVoxEstimate/Temp/R2mask_${sess}.nii.gz \
        ./BestVoxEstimate/Temp/mSD00${sess}.nii.gz

    angstr="${angstr} ./BestVoxEstimate/Temp/mANG00${sess}.nii.gz"
    eccstr="${eccstr} ./BestVoxEstimate/Temp/mECC00${sess}.nii.gz"
    rfsstr="${rfsstr} ./BestVoxEstimate/Temp/mRFS00${sess}.nii.gz"
    sdstr="${sdstr} ./BestVoxEstimate/Temp/mSD00${sess}.nii.gz"
done

echo "Merging "${angstr} 
fslmerge -t ./BestVoxEstimate/Temp/ANGcat.nii.gz ${angstr}
echo "Merging "${eccstr} 
fslmerge -t ./BestVoxEstimate/Temp/ECCcat.nii.gz ${eccstr}
echo "Merging "${rfsstr} 
fslmerge -t ./BestVoxEstimate/Temp/RFScat.nii.gz ${rfsstr}
echo "Merging "${sdstr} 
fslmerge -t ./BestVoxEstimate/Temp/SDcat.nii.gz ${sdstr}

fslmaths ./BestVoxEstimate/Temp/ANGcat.nii.gz -Tmax ./BestVoxEstimate/ANG_bestfit.nii.gz
fslmaths ./BestVoxEstimate/Temp/ECCcat.nii.gz -Tmax ./BestVoxEstimate/ECC_bestfit.nii.gz
fslmaths ./BestVoxEstimate/Temp/RFScat.nii.gz -Tmax ./BestVoxEstimate/RFS_bestfit.nii.gz
fslmaths ./BestVoxEstimate/Temp/SDcat.nii.gz -Tmax ./BestVoxEstimate/SD_bestfit.nii.gz
