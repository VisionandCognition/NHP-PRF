function AverageResults_LISA(MONKEY,R2thr,ReSaveAll)
% Averages results for pRF fits per session
% AverageResults(R2thr,ReSaveAll)
% <MONKEY> specifies which animal
% <R2thr> is the R2 goodness of fit threshold for voxels to be included
%   when done in non-compiled version this can be a vector in which case it
%   will be done for all values, for compiled it should be a scalar
% <ReSaveAll> if true loads individual session files and saves a combined
% file (this also happens when a combined file doesn't exist yet)

%% place definitions & init
if nargin<1
    fprintf('No monkey specified, going with Danny\n');
    MONKEY = 'danny';
    fprintf('No R2 threshold defined. Picking 5\n');
    R2thr = 5;
    ReSaveAll=false;
elseif nargin < 2
    fprintf('No R2 threshold defined. Picking 5\n');
    R2thr = 5;
    ReSaveAll=false;
elseif nargin < 3
    ReSaveAll=false;
end

homefld = pwd;
ResFld =  ['~/Documents/MRI_ANALYSIS/pRF-NHP/Results/' MONKEY '/'];
% ResFld =  ['~/Documents/MRI_ANALYSIS/pRF-NHP/Results/LISA/' MONKEY '/'];

cd(ResFld);
fld = dir('2*');

% convert inputs to correct format
% necessary for compiled
if ischar(R2thr)
    R2thr=eval(R2thr);
end
if ischar(ReSaveAll)
    ReSaveAll=eval(ReSaveAll);
end

if ReSaveAll
    fprintf('ReSaveAll status is TRUE\n');
else
    fprintf('ReSaveAll status is FALSE\n');
end

% Add nifti reading toolbox
tool_basepath = '~/Dropbox/MATLAB_NONGIT/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

%% load the data
R2=[];ANG=[];ECC=[];RFS=[];
if exist('AllFitResults_SD.mat','file') && ~ReSaveAll
    fprintf('Loading existing session data...\n');
    load('AllFitResults_SD');
    %load('AllFitResults_us_motreg');
    % do not resave
    SaveCombined = false;
else
    fprintf('Loading individual session data...\n');
    for s=1:length(fld)
        cd(fullfile(fld(s).folder, fld(s).name));
        matfile = ls('*mat');
        fprintf(['Loading ' matfile '\n']);
        R(s) = load(matfile(1:end-4));
    end
    SaveCombined = true;
end

% reorder
for s=1:length(R)
    R2=cat(4,R2,R(s).result.R2);
    ANG=cat(4,ANG,R(s).result.ang);
    ECC=cat(4,ECC,R(s).result.ecc/10);
    RFS=cat(4,RFS,R(s).result.rfsize/10);
end


%% Loop over thresholds
for tv=1:length(R2thr)
    fprintf(['R2thr is ' num2str(R2thr(tv)) '\n']);
    
    %% threshold based on R2
    thrANG = ANG; thrANG(R2<R2thr(tv)) = NaN;
    % to rad - unwrap - mean - to deg
    thrANG = thrANG*(pi/180); %rad
    thrANG = unwrap(thrANG,[],4); %unwrap
    mtANG = nanmean(thrANG,4); %average
    sdANG = nanstd(thrANG,[],4); %sd
    mtANG = mtANG*(180/pi); %deg
    mtANG(mtANG<0) = mtANG(mtANG<0)+360;
    mtANG(mtANG>360) = mtANG(mtANG>360)-360;
    
    mtR2 = nanmean(R2,4); %average
    sdR2 = nanstd(R2,[],4); %sd


    thrECC = ECC; thrECC(R2<R2thr(tv)) = NaN;
    mtECC = nanmean(thrECC,4);
    sdECC = nanstd(thrECC,[],4);
    
    thrRFS = RFS; thrRFS(R2<R2thr(tv)) = NaN;
    mtRFS = nanmean(thrRFS,4);
    sdRFS = nanstd(thrRFS,[],4);
    
    pVox = mean(R2>R2thr(tv),4);
    
    %% replace NaN with -99
    mtANG(isnan(mtANG)) = -99;
    mtECC(isnan(mtECC)) = -99;
    mtRFS(isnan(mtRFS)) = -99;
    mtR2(isnan(mtR2)) = -99;
    sdANG(isnan(sdANG)) = -99;
    sdECC(isnan(sdECC)) = -99;
    sdRFS(isnan(sdRFS)) = -99;
    sdR2(isnan(sdR2)) = -99;
    
    %% create nifti-files
    cd(ResFld);
    warning off
    mkdir('AveragedResults_SD');cd('AveragedResults_SD');
    mkdir(['Thr_' num2str(R2thr(tv))]);cd(['Thr_' num2str(R2thr(tv))]);
    
    nii = make_nii(mtANG,[1 1 1],[],16,'pRF fit: Angle (deg)');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['MeanAngle_th' num2str(R2thr(tv)) '.nii']);
    gzip(['MeanAngle_th' num2str(R2thr(tv)) '.nii']);
    delete(['MeanAngle_th' num2str(R2thr(tv)) '.nii']);
    
    nii = make_nii(sdANG,[1 1 1],[],16,'pRF fit: Angle (deg)');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['StdAngle_th' num2str(R2thr(tv)) '.nii']);
    gzip(['StdAngle_th' num2str(R2thr(tv)) '.nii']);
    delete(['StdAngle_th' num2str(R2thr(tv)) '.nii']);
    
    nii = make_nii(mtECC,[1 1 1],[],16,'pRF fit: Eccentricity (deg)');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['MeanEccentricity_th' num2str(R2thr(tv)) '.nii']);
    gzip(['MeanEccentricity_th' num2str(R2thr(tv)) '.nii']);
    delete(['MeanEccentricity_th' num2str(R2thr(tv)) '.nii']);
    
    nii = make_nii(sdECC,[1 1 1],[],16,'pRF fit: Eccentricity (deg)');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['StdEccentricity_th' num2str(R2thr(tv)) '.nii']);
    gzip(['StdEccentricity_th' num2str(R2thr(tv)) '.nii']);
    delete(['StdEccentricity_th' num2str(R2thr(tv)) '.nii']);
    
    nii = make_nii(mtRFS,[1 1 1],[],16,'pRF fit: pRF size (deg)');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['MeanRFS_th' num2str(R2thr(tv)) '.nii']);
    gzip(['MeanRFS_th' num2str(R2thr(tv)) '.nii']);
    delete(['MeanRFS_th' num2str(R2thr(tv)) '.nii']);
    
    nii = make_nii(sdRFS,[1 1 1],[],16,'pRF fit: pRF size (deg)');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['StdRFS_th' num2str(R2thr(tv)) '.nii']);
    gzip(['StdRFS_th' num2str(R2thr(tv)) '.nii']);
    delete(['StdRFS_th' num2str(R2thr(tv)) '.nii']);
    
    nii = make_nii(mtR2,[1 1 1],[],16,'pRF fit: R2');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['MeanR2_th' num2str(R2thr(tv)) '.nii']);
    gzip(['MeanR2_th' num2str(R2thr(tv)) '.nii']);
    delete(['MeanR2_th' num2str(R2thr(tv)) '.nii']);
    
    nii = make_nii(sdR2,[1 1 1],[],16,'pRF fit: R2');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['StdR2_th' num2str(R2thr(tv)) '.nii']);
    gzip(['StdR2_th' num2str(R2thr(tv)) '.nii']);
    delete(['StdR2_th' num2str(R2thr(tv)) '.nii']);

    nii = make_nii(pVox,[1 1 1],[],16,'Proportion of sessions in voxel avg');
    nii.hdr.hist = R(1).result.ref_hdr.hist;
    save_nii(nii, ['pVox_th' num2str(R2thr(tv)) '.nii']);
    gzip(['pVox_th' num2str(R2thr(tv)) '.nii']);
    delete(['pVox_th' num2str(R2thr(tv)) '.nii']);
    
    cd ..; cd ..;
end

%% save
cd(ResFld);
if SaveCombined
    fprintf('Saving the combined results file...\n');
    %save('AllFitResults_us_motreg','R','-v7.3');
    save('AllFitResults_SD','R','-v7.3');
end
cd(homefld);

%% Back to start folder
fprintf('Done!\n');