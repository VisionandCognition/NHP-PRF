function ck_GetMRI_pRF(monkeys,models);

%% Add nifti reading toolbox ==============================================
tool_basepath = '~/Dropbox/MATLAB_NONGIT/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

%% Set fMRI paths =========================================================
base_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP-MRI/Projects/pRF/FitResults';
roi_path = ...
     '/Users/chris/Dropbox/CURRENT_PROJECTS/PRF/fMRI_ePhys/Data/fmri/NonReg/ROI';

%% collect the fit results ================================================
for m = 1:length(monkeys)
    frintf(['Processing monkey: ' monkeys{m} '\n']);
    cMonkey = monkeys{m}; sMonkey(1)=upper(cMonkey(1));
    for mm = 1:length(models)
        fprintf(['Model: ' models{mm} '\n']);
        
        RES = load(fullfile(base_path,monkeys{m},models{mm}),'result');
        
        R(m).monkey = monkeys{m};
        R(m).model(mm).fit = models{mm};
        R(m).model(mm).mode = 'MRI';
        R(m).model(mm).hrf = RES.result.options.hrf;
        R(m).model(mm).voldim = size(RES.result.R2,4);
                
        for i = 1:size(RES.result.R2,4)
            R(m).model(mm).ang(i) = reshape(RES.result.ang(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).ecc(i) = reshape(RES.result.ecc(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).rfs(i) = reshape(RES.result.rfs(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).gain(i) = reshape(RES.result.gain(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).R2(i) = reshape(RES.result.R2(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            R(m).model(mm).fwhm(i) = reshape(RES.result.fwhm(:,:,:,i),...
                [numel(RES.result.ang(:,:,:,i)),1]);
            if monkeys{m}(1:3) == 'dog'
                R(m).model(mm).sdratio(i) = reshape(RES.result.sdratio(:,:,:,i),...
                    [numel(RES.result.ang(:,:,:,i)),1]);
                R(m).model(mm).normamp(i) = reshape(RES.result.normamp(:,:,:,i),...
                    [numel(RES.result.ang(:,:,:,i)),1]);
            else
                R(m).model(mm).sdratio(i) = [];
                R(m).model(mm).normamp(i) = [];
            end
            if monkeys{m}(1:3) == 'css'
                R(m).model(mm).expt(i) = reshape(RES.result.expt(:,:,:,i),...
                    [numel(RES.result.ang(:,:,:,i)),1]);
            else
                R(m).model(mm).expt(i) = [];
            end
        end
    end
    clear RES
end

%% Get ROI info ===========================================================
for m = 1:length(monkeys)
    cMonkey = monkeys{m}; sMonkey(1)=upper(cMonkey(1));
    frintf(['Processing ROIs for monkey: ' monkeys{m} '\n']);
    
    ROIs = {'V1', 'V2_merged','V3_merged','V3A','V4_merged','VIP','5_merged',...
        '7_merged','LIP_merged','MST','MT','TEO','TPO',...
        [cMonkey '_LH_V1_electrodes'],[cMonkey '_LH_V4_electrodes']);
    
    for r=1:length(ROIs)
        R(m).ROI(r).name = ROIs{r};
        roiidx = load_nii(...
            fullfile(MRI_path,monkey{m},ROIs{r}));
        
        R(m).ROI(r).idx = reshape(roiidx, [numel(roiidx),1]);
    end

end

