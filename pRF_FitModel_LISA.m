function pRF_FitModel_LISA(Monkey,Session)
% This is the script that runs the pRF fit
% it should be compiled and called from run_CompiledMatlab_LISA.sh
% This means it cannot use 'addpath', instead toolboxes or required
% function should be added when it is compiled usin -a /toolbox etc.

% The uncompiled version of this function should be in
% $HOME/PRF/Code

% this compiled function will be running from:
% "$TMPDIR"/PRF

% Data will be in
% "$TMPDIR"/PRF

% MAKE SURE TO:
% save results in "$TMPDIR"/PRF/Results/%MONKEY/$SESS


% fits the prf model to voxels
% Monkey: string, no caps
% Session: string YYYYMMDD

% These are fixed for this configuration ===
TR=2.5;
doUpsample=true;
doExtraRegression=true;
fitOnlyPosterior=false;
numWorkers =[];
% ==========================================

% Notification of the fact that we're starting
disp(['Starting script for job ' Monkey ',Ses-' Session])

if ~exist(['"$TMPDIR"/PRF/Results/' Monkey '/' Session],'dir')
    mkdir(['"$TMPDIR"/PRF/Results/' Monkey '/' Session])
end

% Link to the brain mask
BrainMask_file = '"$TMPDIR"/PRF/T1_to_func_brainmask_zcrop.nii';

% make outputfolder
result_folder = ['"$TMPDIR"/PRF/Results/' Monkey '/' Session];
if ~exist(result_folder,'dir')
    mkdir(result_folder);
end

%% MODEL PRFs /SESSION ====================================================
% Do the pRF model fit on a session/day basis. Concatenate runs.
% Number of parallel processes might need to be restricted beause of memory
% limitations. Set numWorkers at the start of this function.
% Modeling everything together may be overkill (size wise).

% get the brain mask ----
fprintf('Unpacking BrainMask');
mask_nii=load_nii(BrainMask_file);%load_nii(uz_nii{1});
if fitOnlyPosterior
    mask_nii.img(:,50:end,:)=0; % mask out anterior part (speed up fits)
end
fprintf(' ...done\n');

% run the model-fits
fprintf(['=== Fitting pRF model for ses-' Session ' ===\n']);
% load data -----
fprintf('Loading data...\n');
load(fullfile('"$TMPDIR"/PRF', ['ses-' Session])); % s_run
% concatenate -----
stimulus={};fmri_data={};
fprintf('Concatenating stimuli and volumes...\n');
for r=1:length(s_run)
    stimulus{r}=[]; fmri_data{r}=[];
    for voln = 1:size(s_run(r).vol,2)
        stimulus{r} = cat(3, stimulus{r}, s_run(r).stim{voln}); %#ok<*AGROW>
        fmri_data{r} = cat(4, fmri_data{r}, s_run(r).vol{voln});
    end
    if doExtraRegression
        extraregr{r} = cat(2,s_run(r).motion.estimates,s_run(r).rew);
        % 12 motion correction parameters + reward events
    end
end
    
% fit pRF -----
% get indices to mask voxels > 0
options.vxs = find(mask_nii.img>0);
% add regressors when wanted
if doExtraRegression
    options.wantglmdenoise = extraregr;
end
    
% start a parallel pool of workers
if ~isempty(numWorkers) && s==1
    parpool(numWorkers)
else
    % don't predefine the number of workers
    % let it take the max available when running
end

% run analyzePRF tool
if doUpsample % tr = TR/2
    result = analyzePRF(stimulus,fmri_data,TR/2,options);
else
    result = analyzePRF(stimulus,fmri_data,TR,options);
end

% save the result ----
fprintf('Saving the result: ');
save(fullfile(result_folder,['pRF_Sess-' Session]),...
    'result','-v7.3');
% also save as nifti files
% angle ---
fprintf('Angles ');
nii = make_nii(result.ang,[1 1 1],[],[],...
    'pRF fit: Angles (deg)');
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_ang.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_ang.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_ang.nii']));
% ecc ---
fprintf('Ecc ');
nii = make_nii(result.ecc,[1 1 1],[],[],...
    'pRF fit: Eccentricity (pix)');
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_ecc.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_ecc.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_ecc.nii']));
% size ---
fprintf('Size ');
nii = make_nii(result.rfsize,[1 1 1],[],[],...
    'pRF fit: RF size (pix)');
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_rfsize.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_rfsize.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_rfsize.nii']));
% R^2 Goodness of fit ---
fprintf('R2 ');
nii = make_nii(result.R2,[1 1 1],[],[],...
    'pRF fit: R2 Goodnes off fit');
save_nii(nii, fullfile(result_folder, ['Sess-' Session '_R2.nii']));
gzip(fullfile(result_folder, ['Sess-' Session '_R2.nii']));
delete(fullfile(result_folder, ['Sess-' Session '_R2.nii']));

fprintf('>> Done!\n');