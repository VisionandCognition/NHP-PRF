function pRF_FitModel_LISA(Monkey,Session,Slices,numWorkers)
% This is the script that runs the pRF fit
% it should be compiled and called from run_CompiledMatlab_LISA.sh
% This means it cannot use 'addpath', instead toolboxes or required
% function should be added when it is compiled usin -a /toolbox etc.

% complie on LISA with:
% module load matlab
% echo "mcc -m $HOME/PRF/Code/pRF_FitModel_LISA.m -a $HOME/PRF/Code/analyzePRF -a $HOME/PRF/Code/NIfTI" | matlab -nodisplay
% module unload matlab

% The uncompiled version of this function should be in
% $HOME/PRF/Code

% this compiled function will be running from:
% "$TMPDIR"/PRF

% Data will be in
% "$TMPDIR"/PRF

% MAKE SURE TO:
% save results in "$TMPDIR"/PRF/Results/%MONKEY/$SESS/$SLICES

% fits the prf model to voxels
% Monkey: string, no caps
% Session: string YYYYMMDD

% These are fixed for this configuration ===
TR=2.5;
if strcmp(Monkey, 'eddy') && strcmp(Session,'20160721')
    TR=3;
end

doUpsample=true;
doExtraRegression=true;
mlroot = pwd; % $TMPDIR/PRF
% ==========================================

if ischar(Slices)
    slices = eval(Slices);
else
    slices=Slices;
end

if ischar(numWorkers)
    numWorkers = eval(numWorkers);
end


SliceLabel = [num2str(slices(1),'%02d') '-' num2str(slices(end),'%02d')];

% Notification of the fact that we're starting
disp(['Starting script for job ' Monkey ', Ses-' Session ', Slices ' Slices])

% Link to the brain mask
if strcmp(Monkey, 'danny')
    BrainMask_file = [mlroot '/T1_to_func_brainmask_zcrop.nii'];
elseif strcmp(Monkey, 'eddy')
    BrainMask_file = [mlroot '/HiRes_to_T1_mean.nii_shadowreg_Eddy_brainmask.nii'];
end

% make outputfolder
result_folder = [mlroot '/Results/' Monkey '/' Session '/Slices_' SliceLabel];
if ~exist(result_folder,'dir')
    mkdir(result_folder);
end
fprintf(['Saving results in: ' result_folder]);

%% MODEL PRFs /SESSION ====================================================
% Do the pRF model fit on a session/day basis. Concatenate runs.
% Number of parallel processes might need to be restricted beause of memory
% limitations. Set numWorkers at the start of this function.
% Modeling everything together may be overkill (size wise).

% get the brain mask ----
fprintf('\nUnpacking BrainMask');
mask_nii=load_nii(BrainMask_file);%load_nii(uz_nii{1});
fprintf(' ...done\n');

% run the model-fits
fprintf(['=== Fitting pRF model for ses-' Session ' ===\n']);
% load data -----
fprintf('Loading data...\n');
load(fullfile(mlroot, ['ses-' Session])); % s_run
% concatenate -----
stimulus={};fmri_data={};
fprintf('Concatenating stimuli and volumes...\n');
for r=1:length(s_run)
    stimulus{r}=[]; fmri_data{r}=[];
    for voln = 1:size(s_run(r).vol,2)
        stimulus{r} = cat(3, stimulus{r}, s_run(r).stim{voln}); %#ok<*AGROW>
        %fmri_data{r} = cat(4, fmri_data{r}, s_run(r).vol{voln});
        fmri_data{r} = cat(4, fmri_data{r}, s_run(r).vol{voln}(:,:,slices));
    end
    if doExtraRegression
        extraregr{r} = cat(2,s_run(r).motion.estimates,s_run(r).rew);
        % 12 motion correction parameters + reward events
    end
end
    
% fit pRF -----
% get indices to mask voxels > 0
options.vxs = find(mask_nii.img(:,:,slices)>0);
options.display = 'final';

% add regressors when wanted
if doExtraRegression
    options.wantglmdenoise = extraregr;
end
    
% start a parallel pool of workers
if ~isempty(numWorkers)
    parpool(numWorkers)
else
    % if numWorkers = []
    % don't predefine the number of workers
    % let it take the max available when running
end

% run analyzePRF tool
if doUpsample % tr = TR/2
    result = analyzePRF(stimulus,fmri_data,TR/2,options);
else
    result = analyzePRF(stimulus,fmri_data,TR,options);
end
result.hdr_ref = s_run(1).hdr;

% save the result ----
fprintf('\nSaving the result: ');
save(fullfile(result_folder,['pRF_Sess-' Session]),...
    'result','-v7.3');
% also save as nifti files
% angle ---
fprintf('Angles ');
nii = make_nii(result.ang,[1 1 1],[],32,...
    'pRF fit: Angles (deg)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_ang.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_ang.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_ang.nii']));
% ecc ---
fprintf('Ecc ');
nii = make_nii(result.ecc,[1 1 1],[],32,...
    'pRF fit: Eccentricity (pix)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_ecc.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_ecc.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_ecc.nii']));
% size ---
fprintf('Size ');
nii = make_nii(result.rfsize,[1 1 1],[],32,...
    'pRF fit: RF size (pix)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_rfsize.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_rfsize.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_rfsize.nii']));
% R^2 Goodness of fit ---
fprintf('R2 ');
nii = make_nii(result.R2,[1 1 1],[],32,...
    'pRF fit: R2 Goodnes off fit');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_R2.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_R2.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_R2.nii']));

cd ..
fprintf('>> Done!\n');