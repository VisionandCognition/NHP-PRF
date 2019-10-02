#!/bin/bash

# register template epi (space of results) to the anatomical in NMT single subjects
# for this anatomical space we sill also have the surfaces 
# (but these require additional registration because of the header manipulation 
# that is required to get Freesurfer working)


SUBJ=$1  # give    monkey name as argument

if [ $SUBJ == 'danny' ]; then
    t1_brain=/NHP_MRI/Template/NMT/single_subject_scans/Danny/NMT_Danny_process/Danny_brain.nii.gz
    wm=/NHP_MRI/Template/NMT/single_subject_scans/Danny/NMT_Danny_process/Danny_segmentation_WM.nii.gz
    epi=../Results/Reference/$SUBJ/input_files/Func/sub-danny_ses-20180117_task-prf_run-01_frame-10_bold_res-1x1x1_reference_zcrop.nii.gz
    epi_brainmask=../Results/Reference/$SUBJ/input_files/Func/T1_to_func_brainmask_zcrop.nii
elif [ $SUBJ == 'eddy' ]; then
    t1_brain=/NHP_MRI/Template/NMT/single_subject_scans/Eddy/NMT_Eddy_process/Eddy_brain.nii.gz
    wm=/NHP_MRI/Template/NMT/single_subject_scans/Eddy/NMT_Eddy_process/Eddy_segmentation_WM.nii.gz
    epi=../Results/Reference/$SUBJ/input_files/Func/ref_func_undist_inData_al_fnirt.nii.gz
    epi_brainmask=../Results/Reference/$SUBJ/input_files/Func/HiRes_to_T1_mean_shadowreg_Eddy_brainmask.nii.gz
fi

# mask out epi-brain
fslmaths $epi -mas $epi_brainmask ../Results/Reference/$SUBJ/input_files/Func/epi_brain.nii.gz
epi_brain=../Results/Reference/$SUBJ/input_files/Func/epi_brain.nii.gz

# regster epi to t1
if [ $SUBJ == 'danny' ]; then
    flirt -ref ${t1_brain} -wmseg ${wm} -in ${epi-brain} \
        -out ../Results/Reference/$SUBJ/output_files/epi2anat.nii.gz -omat ../Results/Reference/$SUBJ/output_files/epi2anat.mat -pedir -2   
elif [ $SUBJ == 'eddy' ]; then
    flirt -dof 6 \
        -ref ${t1_brain} -wmseg ${wm} -in ${epi-brain} \
        -out ../Results/Reference/$SUBJ/output_files/epi2anat.nii.gz -omat ../Results/Reference/$SUBJ/output_files/epi2anat.mat -pedir -2   
fi

# use the reg matrix to warp the results
mkdir -p ../Results/$SUBJ/Reg2Anat/BestVoxEstimate

flirt -ref ${t1_brain} \
    -in ../Results/$SUBJ/BestVoxEstimate/ANG_bestfit.nii.gz \
    -applyxfm -init ../Results/Reference/$SUBJ/output_files/epi2anat.mat \
    -interp nearestneighbour  \
    -out ../Results/$SUBJ/Reg2Anat/BestVoxEstimate/ANG.nii.gz &

flirt -ref ${t1_brain} \
    -in ../Results/$SUBJ/BestVoxEstimate/ECC_bestfit.nii.gz \
    -applyxfm -init ../Results/Reference/$SUBJ/output_files/epi2anat.mat \
    -interp nearestneighbour  \
    -out ../Results/$SUBJ/Reg2Anat/BestVoxEstimate/ECC.nii.gz &

flirt -ref ${t1_brain} \
    -in ../Results/$SUBJ/BestVoxEstimate/RFS_bestfit.nii.gz \
    -applyxfm -init ../Results/Reference/$SUBJ/output_files/epi2anat.mat \
    -interp nearestneighbour  \
    -out ../Results/$SUBJ/Reg2Anat/BestVoxEstimate/RFS.nii.gz &

flirt -ref ${t1_brain} \
    -in ../Results/$SUBJ/BestVoxEstimate/SD_bestfit.nii.gz \
    -applyxfm -init ../Results/Reference/$SUBJ/output_files/epi2anat.mat \
    -interp nearestneighbour  \
    -out ../Results/$SUBJ/Reg2Anat/BestVoxEstimate/SD.nii.gz &

flirt -ref ${t1_brain} \
    -in ../Results/$SUBJ/BestVoxEstimate/R2_max.nii.gz \
    -applyxfm -init ../Results/Reference/$SUBJ/output_files/epi2anat.mat \
    -interp nearestneighbour  \
    -out ../Results/$SUBJ/Reg2Anat/BestVoxEstimate/R2.nii.gz &

fslmaths ../Results/$SUBJ/BestVoxEstimate/brainmask.nii.gz -thr 1 -bin ../Results/$SUBJ/BestVoxEstimate/brainmask_tight.nii.gz

flirt -ref ${t1_brain} \
    -in ../Results/$SUBJ/BestVoxEstimate/brainmask_tight.nii.gz \
    -applyxfm -init ../Results/Reference/$SUBJ/output_files/epi2anat.mat \
    -interp nearestneighbour  \
    -out ../Results/$SUBJ/Reg2Anat/BestVoxEstimate/brainmask.nii.gz &