function pRF_FitModel_LISA(Monkey,Session)
% fits the prf model to voxels
% Monkey: string, no caps
% Session: cell array with YYYYMMDD
% NB inputs are generated in other script
TR=2.5;



doUpsample=true;

doExtraRegression=true;

fitOnlyPosterior=false;


%% Add toolbox paths <<< CAN'T DO THIS IN COMPILED CODE


tool_basepath = '/Users/chris/Documents/MATLAB/TOOLBOX';
BIDS_basepath = '/NHP_MRI/NHP-BIDS/';

addpath(genpath('/media/DOCUMENTS/DOCUMENTS/MRI_ANALYSIS/analyzePRF'));


% Add nifti reading toolbox
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));
% Add Kendrick Kay's pRF analysis toolbox
addpath(genpath(fullfile(tool_basepath, 'analyzePRF')));

% Link to the brain mask
if strcmp(MONKEY, 'danny')
    BrainMask_file = fullfile(BIDS_basepath, 'manual-masks','sub-danny',...
        'ses-20180117','func','T1_to_func_brainmask_zcrop.nii');
else
    error('Unknown monkey name or no mask available')
end

% create a folder to save outputs in
if doUpsample
    out_folder = ['pRF_sub-' MONKEY '_us'];
else
    out_folder = ['pRF_sub-' MONKEY]; %#ok<*UNRCH>
end
warning off %#ok<*WNOFF>
mkdir(out_folder);
warning on %#ok<*WNON>





%% MODEL PRFs /SESSION ====================================================
% Do the pRF model fit on a session/day basis. Concatenate runs.
% Number of parallel processes might need to be restricted beause of memory
% limitations. Set numWorkers at the start of this function.
% Modeling everything together may be overkill (size wise).

% outputfolder ------
if doUpsample
    if doExtraRegression
        result_folder = ['FitResult_motregr_sub-' MONKEY];
    else
        result_folder = ['FitResult_us_sub-' MONKEY];
    end
else
    result_folder = ['FitResult_sub-' MONKEY];
end
warning off; mkdir(result_folder); warning on;

% get the brain mask ----
fprintf('Unpacking BrainMask');
mask_nii=load_nii(BrainMask_file);%load_nii(uz_nii{1});
if fitOnlyPosterior
    mask_nii.img(:,50:end,:)=0; % mask out anterior part (speed up fits)
end
fprintf(' ...done\n');

% run the model-fits
for s=1:length(sessions)
    fprintf(['=== Fitting pRF model for ses-' sessions{s} ' ===\n']);
    
    % load data -----
    fprintf('Loading data...\n');
    load(fullfile(out_folder, ['ses-' sessions{s}])); % s_run
    
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
    
    if doUpsample % tr = TR/2
        result = analyzePRF(stimulus,fmri_data,TR/2,options);
    else
        result = analyzePRF(stimulus,fmri_data,TR,options);
    end
    
    % save the result ----
    fprintf('Saving the result: ');
    save(fullfile(result_folder,['pRF_Sess-' sessions{s}]),...
        'result','-v7.3');
    % also save as nifti files
    % angle ---
    fprintf('Angles ');
    nii = make_nii(result.ang,[1 1 1],[],[],...
        'pRF fit: Angles (deg)');
    save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s} '_ang.nii']));
    gzip(fullfile(result_folder, ['Sess-' sessions{s} '_ang.nii']));
    delete(fullfile(result_folder, ['Sess-' sessions{s} '_ang.nii']));
    % ecc ---
    fprintf('Ecc ');
    nii = make_nii(result.ecc,[1 1 1],[],[],...
        'pRF fit: Eccentricity (pix)');
    save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s} '_ecc.nii']));
    gzip(fullfile(result_folder, ['Sess-' sessions{s} '_ecc.nii']));
    delete(fullfile(result_folder, ['Sess-' sessions{s} '_ecc.nii']));
    % size ---
    fprintf('Size ');
    nii = make_nii(result.rfsize,[1 1 1],[],[],...
        'pRF fit: RF size (pix)');
    save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s} '_rfsize.nii']));
    gzip(fullfile(result_folder, ['Sess-' sessions{s} '_rfsize.nii']));
    delete(fullfile(result_folder, ['Sess-' sessions{s} '_rfsize.nii']));
    % R^2 Goodness of fit ---
    fprintf('R2 ');
    nii = make_nii(result.R2,[1 1 1],[],[],...
        'pRF fit: R2 Goodnes off fit');
    save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s} '_R2.nii']));
    gzip(fullfile(result_folder, ['Sess-' sessions{s} '_R2.nii']));
    delete(fullfile(result_folder, ['Sess-' sessions{s} '_R2.nii']));
    
    fprintf('>> Done!\n');
end




















%% ENNY'S STUFF >>>>>>>>


%% set path to your m-files
disp(['Starting script for job nr ' num2str(job_nr1) '_' num2str(job_nr2) '_' num2str(job_nr3) '_' num2str(job_nr4)])

% example 1: add current path
curr_path = pwd;
addpath(curr_path)  % adds path only for this session

if ~exist('"$TMPDIR"/PRF/TMPResults','dir')
    mkdir('"$TMPDIR"/TMPResults')
end
addpath(genpath('"$TMPDIR"/TMPData')) % Set path
addpath(genpath('"$TMPDIR"/TMPResults')) % Set path
addpath(genpath('"$TMPDIR"/MVPA_Scripts')) %Set path

%% copy your data from the server to the local machine
%% run your script
% if you copy files, make sure to surround your script with a try/catch
% loop, so that you can still copy the resulting files back
try
    %     load('/home/enny/Documents/ADNI_NEW/TemporarySaveDir.mat');
    %     loadfile = [tempsave 'TemporarySave2']
    %     load(loadfile)
    if ischar(job_nr2)
        job_nr2 = str2num(job_nr2);
    end
    if ischar(job_nr3)
        job_nr3 = str2num(job_nr3);
    end
    if ischar(job_nr4)
        job_nr4 = str2num(job_nr4);
    end
    display(['Running script for job nr ' num2str(job_nr1) '_' num2str(job_nr2) '_' num2str(job_nr3) '_' num2str(job_nr4)])
    % call your script here, e.g.
    %         RunModel(job_nr, loadfile, X, T, nfolds)
    SVM_StimDecoder_AcrossTimeWindows(job_nr1,job_nr2,job_nr3,job_nr4)
    %%pRF_FitModel_core_LISA(Monkey,Sessions,doUpsample,doExtraRegression,fitOnlyPosterior)
catch my_error
    % display error and stack, and continue with the rest of the script
    display_error_info(my_error)
    display(['Script for job nr ' num2str(job_nr1) '_' num2str(job_nr2) '_' num2str(job_nr3) '_' num2str(job_nr4) ' aborted'])
end


%% copy you data back and clean up

% in case you are using a local directory, copy all data back to the net,
% and dont forget to delete the data on the remote machine afterwards

display(['Finished script for job nr ' num2str(job_nr1) '_' num2str(job_nr2) '_' num2str(job_nr3) '_' num2str(job_nr4)])

end

function display_error_info(ME)
ME
for curr_field = fieldnames(ME)'
    display([curr_field{1} ':'])
    display(ME.(curr_field{1}))
    display(' ')
end
for stack_idx = 1:size(ME.stack, 1)
    display(ME.stack(stack_idx))
    display(' ')
end
end