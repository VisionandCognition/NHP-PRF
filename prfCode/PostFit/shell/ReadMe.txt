====================================
Workflow POSTFIT pRF's (shell)
====================================

These post-fitting scripts are located in the <..>/prfCode/PostFit/shell folder

NB! Functionality of these scripts has been replace by Jupyter Notebooks
<..>/prfCode./Notebooks/MRI_pRF_Process.ipynb

====================================

<..>/CrossVal/Batch_Mask_MRI_Result_R2_cv1.sh
	- runs Mask_MRI_Result_R2_cv1.sh in a loop
		- Processses MRI results and masks based on R2 
		- Gets results data from <..>/FitResults/MRI/<subject>/<model>
		- Gets a brainmask file from <..>/FitResults/Reference/Volumes/func/${MONKEY}/sub-${MONKEY,}_ref_func_mask_res-1x1x1.nii.gz

<..>/Generic/Mask_with_ROI.sh
	- Masks result volumes with an ROI
	-  
