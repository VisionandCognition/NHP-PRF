#!/bin/bash
MONKEY=$1
MODEL=$2
ROI=$3

echo Processing monkey ${MONKEY} 
echo ${MODEL} ${ROI} 

M_FLD=../../../../FitResults/MRI/${MONKEY}/${MODEL}/inAnat
ROI_FLD=../../../../FitResults/Reference/Volumes/atlas/${MONKEY}/ROI_manualadjust
home_fld=${pwd}

# do this for all existing thresholded subfolders 
for sf in ${M_FLD}/TH_*; do
	if [ -d ${sf} ]; then
		echo 'Processing '$sf
		
		# mask averaged results with these ROIs
		# NB! Use the filenames from ROI_manualadjust
		# use specified ROI or this list

		if [ ${ROI} == 'allrois']; then 
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
	    else;
	    	declare -a arr=(
	    		${ROI}	
	    	)
	    fi

		for roi in "${arr[@]}"; do
			echo 'ROI '${roi}

			mkdir -p ${sf}/mROI/${roi}
			
			fslmaths ${ROI_FLD}/${roi}_roi.nii.gz -sub 1 -mul 99 \
					 ${ROI_FLD}/${roi}_SUB.nii

			# mask with roi file
			for file in ${sf}/*.nii.gz; do
				fslmaths ${file} -mas ${ROI_FLD}/${roi}_roi.nii.gz \
					 ${sf}/mROI/${roi}/$(basename "$file")
				fslmaths ${sf}/mROI/${roi}/$(basename "$file") -add \
						 ${ROI_FLD}/${roi}_SUB.nii.gz \
						 ${sf}/mROI/${roi}/$(basename "$file")
			done
			rm ${ROI_FLD}/${roi}_SUB*
		done
	fi
done
