#!/bin/bash
MONKEY=$1

echo 'Pocessing monkey '${MONKEY} 

if [ ${MONKEY} = 'Danny' ]; then
	M_FLD=../Results/danny/AveragedResults
	ROI_FLD=../Results/Reference/danny/output_files/ROI_adj/1mm/nii
else  # Eddy
	M_FLD=../Results/eddy/AveragedResults
	ROI_FLD=../Results/Reference/eddy/output_files/ROI/1mm/nii
fi

home_fld=${pwd}

# do this for all existing AveragedResults subfolders 
# (based on Thr)
for sf in ${M_FLD}/Thr_*; do
	if [ -d ${sf} ]; then
		echo 'Processing '$sf
		
		# mask averaged results with these ROIs
		declare -a arr=(
    		#'V1'
    		#'V2'
    		#'V3d'
    		#'V3v'
    		#'V3A'
    		#'V4'
    		#'V4t'
    		#'V4v'
    		#'MT'
    		#'MST'
    		#'TPO'
    		#'TEO'
    		#'7a_(Opt-PG)'
    		'5_(PE)'
    		'5_(PEa)'
    		'LIPd'
    		'LIPv'
    		'VIP'	
    		)

		for roi in "${arr[@]}"; do
			echo 'ROI '${roi}

			mkdir -p ${sf}/mROI/
			mkdir -p ${sf}/mROI/${roi}

			
			fslmaths ${ROI_FLD}/${roi}_roi.nii -sub 1 -mul 99 \
					 ${ROI_FLD}/${roi}_SUB.nii

			# mask with roi file
			for file in ${sf}/*.nii.gz; do
				fslmaths ${file} -mas ${ROI_FLD}/${roi}_roi.nii \
					 ${sf}/mROI/${roi}/$(basename "$file")
				fslmaths ${sf}/mROI/${roi}/$(basename "$file") -add \
						 ${ROI_FLD}/${roi}_SUB.nii \
						 ${sf}/mROI/${roi}/$(basename "$file")
			done
			rm ${ROI_FLD}/${roi}_SUB*
		done
	fi
done
