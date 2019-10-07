function ck_BasicProc_AllFits

% this script will do some basic processing like calculate 
% - means from two-way crossvalidated results
% - X,Y coordinates from polar coordinates
% 
% 

fitres_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/';

load(fullfile(fitres_path,'MultiModal','AllFits_cv1'),'R_MRI','R_EPHYS');



rfs
ang
ecc
fwhm
imag
real
maxR2
meanR2



Sess-linhrf_cv1_dhrf_rfsize_2.nii.gz
brainmask.nii.gz
nonbrainmask.nii.gz
pRF_Sess-linhrf_cv1_dhrf.mat
Sess-linhrf_cv1_dhrf_ang_1.nii.gz
Sess-linhrf_cv1_dhrf_ang_2.nii.gz
Sess-linhrf_cv1_dhrf_ecc_1.nii.gz
Sess-linhrf_cv1_dhrf_ecc_2.nii.gz
Sess-linhrf_cv1_dhrf_FWHM_1.nii.gz
Sess-linhrf_cv1_dhrf_FWHM_2.nii.gz
Sess-linhrf_cv1_dhrf_imag_1.nii.gz
Sess-linhrf_cv1_dhrf_imag_2.nii.gz
Sess-linhrf_cv1_dhrf_maxR2.nii.gz
Sess-linhrf_cv1_dhrf_meanR2.nii.gz
Sess-linhrf_cv1_dhrf_minR2.nii.gz
Sess-linhrf_cv1_dhrf_R2_1.nii.gz
Sess-linhrf_cv1_dhrf_R2_2.nii.gz
Sess-linhrf_cv1_dhrf_real_1.nii.gz
Sess-linhrf_cv1_dhrf_real_2.nii.gz
Sess-linhrf_cv1_dhrf_rfsize_1.nii.gz