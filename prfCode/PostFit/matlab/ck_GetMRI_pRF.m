function ck_GetMRI_pRF(monkeys,models,output)

if nargin < 2
    fprintf('ERROR: Not enough arguments specified\n');
    return
end
clc;

%% Add nifti reading toolbox ==============================================
tool_basepath = '~/Dropbox/MATLAB_NONGIT/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

%% Set fMRI paths =========================================================
fitres_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/MRI';
roi_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/Results/MRI';
bids_path = ...
    '/Users/chris/Documents/MRI_ANALYSIS/NHP-BIDS/manual-masks';

%% collect the fit results ================================================
fprintf('Collecting MRI retinotopic maps...\n');

for m = 1:length(monkeys)
    fprintf(['Processing monkey: ' monkeys{m} '\n']);
    cMonkey = monkeys{m}; cMonkey(1)=upper(cMonkey(1));
    
    R(m).monkey = monkeys{m};
    R(m).mode = 'MRI';
    
    for mm = 1:length(models)
        fprintf(['Model: ' models{mm} '\n']);
        
        RES = load(fullfile(fitres_path,monkeys{m},models{mm},...
            ['pRF_Sess-' models{mm}]),'result');
        
        R(m).model(mm).prfmodel = models{mm};
        R(m).model(mm).hrf = RES.result.options.hrf;
        R(m).model(mm).voldim = size(RES.result.R2,4);
                
        for i = 1:size(RES.result.R2,4)
            R(m).model(mm).ang(i,:) = reshape(RES.result.ang(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).ecc(i,:) = reshape(RES.result.ecc(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).rfs(i,:) = reshape(RES.result.rfsize(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).gain(i,:) = reshape(RES.result.gain(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).R2(i,:) = reshape(RES.result.R2(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).fwhm(i,:) = reshape(RES.result.fwhm(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            if models{mm}(1:3) == 'dog'
                R(m).model(mm).sdratio(i,:) = reshape(RES.result.sdratio(:,:,:,i),...
                    [numel(RES.result.ang(:,:,:,i)),1]);
                R(m).model(mm).normamp(i,:) = reshape(RES.result.normamp(:,:,:,i),...
                    [numel(RES.result.ang(:,:,:,i)),1]);
            else
                R(m).model(mm).sdratio = [];
                R(m).model(mm).normamp = [];
            end
            if models{mm}(1:3) == 'css'
                R(m).model(mm).expt(i,:) = reshape(RES.result.expt(:,:,:,i),...
                    [numel(RES.result.ang(:,:,:,i)),1]);
            else
                R(m).model(mm).expt = [];
            end
        end
    end
    clear RES
    
    %% Get ROI info =======================================================
    fprintf(['Processing ROIs for monkey: ' monkeys{m} '\n']);

    ROIs = {'V1', 'V2_merged','V3_merged','V3A','V4_merged','VIP','5_merged',...
        '7_merged','LIP_merged','MST','MT','TEO','TPO',...
        [cMonkey '_LH_V1_electrodes'],[cMonkey '_LH_V4_electrodes']};
    
    BRAINMASK = load_nii(fullfile(bids_path,['sub-' monkeys{m}],'func',...
            ['sub-' monkeys{m} '_ref_func_mask_res-1x1x1.nii']));
    R(m).BRAIN = reshape(BRAINMASK.img, [numel(BRAINMASK.img),1]);
    
    for r=1:length(ROIs)
        fprintf([ROIs{r} '\n']);
        R(m).ROI(r).label = ROIs{r};
        roi_idx = load_nii(fullfile(roi_path, cMonkey,'NonReg',...
            'ROI',[ROIs{r} '.nii']));
        R(m).ROI(r).idx = reshape(roi_idx.img, [numel(roi_idx.img),1]);
    end
end

%% Save the combined results ==============================================
fprintf('Saving the combined MRI result-file\n');
[~,~,~] = mkdir(fullfile(fitres_path,'Combined'));
save(fullfile(fitres_path,'Combined',output),'R')
