%% WHICH DATA =============================================================
Danny_pRF_goodruns;

    
    
%% INITIALIZE =============================================================
% Add nifti reading toolbox
addpath(genpath('/Users/chris/Documents/MATLAB/TOOLBOX/NIfTI'))
% Add Kendrick Kay's pRF analysis toolbox
addpath(genpath('/Users/chris/Documents/MRI_ANALYSIS/analyzePRF'))

%% GET THE FILE-PATHS OF THE IMAGING FILES ================================
% all functional runs that are preprocessed with the BIDS pipeline are 
% resampled to 1x1x1 mm isotropic voxels, reoriented from sphinx,
% motion corrected, (potentially smoothed with 2 mm FWHM), and
% registered to an example functional volume 
% (so they're already in a common space)
% do the analysis in this functional space than we can register to hi-res
% anatomical data and/or the NMT template later
mpath = ['/NHP-MRI/NHP-BIDS/derivatives/featpreproc/highpassed_files/sub-' MONKEY ];
spath = [mpath '/ses-' session '/func'];
rpath = [spath '/' ls([spath '/*run-' runnr '*'])]