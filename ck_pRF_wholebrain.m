function ck_pRF_wholebrain(SessionList, do_Resave, do_FitPRF_perSession)
% collects data, concatenates, downsamples stimulus and resaves
% fits the prf model to voxels

%% WHICH DATA =============================================================
%clear all; clc;
if nargin <3
    fprintf('Not enough arguments specified, will use defaults:\n');
    fprintf('SessionList: Danny_pRF_goodruns\n');
    fprintf('do_Resave = false; doFitPRF_perSession = false\n');
    Danny_pRF_goodruns;
    do_Resave = false;
    do_FitPRF_perSession = false;
else
    eval(SessionList);
end

%% 
TR=2.5;

% This is the 'pure' sweep to volume map
% For the 230 volume version (15 vol blanks)
SwVolMap_new = {    1 , 6:25    ;...
    2 , 41:60   ;...
    3 , 61:80   ;...
    4 , 96:115  ;...
    5 , 116:135 ;...
    6 , 151:170 ;...
    7 , 171:190 ;...
    8 , 206:225 };
% For the 210 volume version (10 vol blanks)
SwVolMap_old = {    1 , 6:25    ;...
    2 , 36:55   ;...
    3 , 56:75   ;...
    4 , 86:105  ;...
    5 , 106:125 ;...
    6 , 136:155 ;...
    7 , 156:175 ;...
    8 , 186:205 };

%% INITIALIZE =============================================================
% Platform specific basepath
if ispc
    tool_basepath = 'D:\CK\code\MATLAB';
    BIDS_basepath = '\\vcnin\NHP_MRI\NHP-BIDS';
else
   tool_basepath = '/Users/chris/Documents/MATLAB/TOOLBOX';
   BIDS_basepath = '/NHP_MRI/NHP-BIDS/';
   addpath(genpath('/media/DOCUMENTS/DOCUMENTS/MRI_ANALYSIS/analyzePRF'));
end
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
out_folder = ['pRF_sub-' MONKEY];
warning off %#ok<*WNOFF>
mkdir(out_folder);
warning on %#ok<*WNON>

%% GET THE FILE-PATHS OF THE IMAGING  & STIM-MASK FILES ===================
% All functional runs that are preprocessed with the BIDS pipeline are
% resampled to 1x1x1 mm isotropic voxels, reoriented from sphinx,
% motion corrected, (potentially smoothed with 2 mm FWHM), and
% registered to an example functional volume
% (so they're already in a common space)
% do the analysis in this functional space than we can register to hi-res
% anatomical data and/or the NMT template later
sessions = unique(DATA(:,1));

monkey_path_nii = fullfile(BIDS_basepath, 'derivatives',...
    'featpreproc','highpassed_files',['sub-' MONKEY]);
monkey_path_stim = fullfile(BIDS_basepath,['sub-' MONKEY]);
for s=1:length(sessions)
    sess_path_nii{s} = fullfile(monkey_path_nii, ['ses-' sessions{s}], 'func'); %#ok<*SAGROW>
    sess_path_stim{s} = fullfile(monkey_path_stim, ['ses-' sessions{s}], 'func');
    runs = unique(DATA(strcmp(DATA(:,1),sessions{s}),2));
    for r=1:length(runs)
        if ispc % the ls command works differently in windows
            a=ls( fullfile(sess_path_nii{s},['*run-' runs{r} '*.nii.gz']));
            run_path_nii{s,r} = fullfile(sess_path_nii{s},a(1:end-3));
            b = ls( fullfile(sess_path_stim{s}, ...
                ['*run-' runs{r} '*model*']));
            run_path_stim{s,r}= fullfile(sess_path_stim{s},b,'StimMask.mat');
        else
            a = ls( fullfile(sess_path_nii{s},['*run-' runs{r} '*.nii.gz']));
            run_path_nii{s,r} = a(1:end-3);
            run_path_stim{s,r} = ls( fullfile(sess_path_stim{s}, ...
                ['*run-' runs{r} '*model*'],'StimMask.mat'));
            run_path_nii{s,r} = run_path_nii{s,r}(1:end-1);
            run_path_stim{s,r} = run_path_stim{s,r}(1:end-1);
        end
        sweepinc{s,r} = DATA( ...
            (strcmp(DATA(:,1),sessions{s}) & strcmp(DATA(:,2),runs{r})),3);
    end
end

%% LOAD & RE-SAVE STIMULUS MASKS & NIFTI ==================================
if do_Resave
    for s=1:size(run_path_stim,1)
        fprintf(['Processing session ' sessions{s} '\n']);
        rps = [];
        for i=1:size(run_path_stim,2) 
            if ~isempty(run_path_stim{s,i}) 
                rps=[rps i]; 
            end 
        end
        for r=rps 
            % stimulus mask -----
            load(run_path_stim{s,r}(1:end-4));
            % loads variable called stimulus (x,y,t) in volumes
            sinc = cell2mat(sweepinc{s,r});
            if size(stimulus,3) == 230
                SwVolMap = SwVolMap_new;
            elseif size(stimulus,3) == 210
                SwVolMap = SwVolMap_old;
            else
                error('weird number of stimulus frames');
            end
            firstvol = SwVolMap{min(sinc),2}(1) - 5;
            lastvol = SwVolMap{max(sinc),2}(end) + 5;
            vinc=firstvol:lastvol;
            % volumes ------
            fprintf('Unpacking nii.gz');
            %uz_nii=gunzip(run_path_nii{s,r});
            temp_nii=load_nii(run_path_nii{s,r});%load_nii(uz_nii{1});
            %delete(uz_nii{1});
            fprintf(' ...done\n');
            % save the session-based stims & vols -----
            for v=1:length(vinc)
                % resample image (160x160 pix gives 10 pix/deg)
                s_run(r).stim{v} = imresize(stimulus(:,:,vinc(v)),[160 160]);
                s_run(r).vol{v} = temp_nii.img(:,:,:,vinc(v));
            end
            clear stimulus temp_nii
        end
        fprintf(['Saving ses-' sessions{s} '\n']);
        save(fullfile(out_folder, ['ses-' sessions{s}]),'s_run','-v7.3');
        clear s_run
    end
else
    fprintf('Not doing the resave (assuming this has already been done)\n');
end

%% MODEL PRFs /SESSION ====================================================
% Do the pRF model fit on a session/day basis. Concatenate runs.
% Modeling everything together may be overkill (size wise)
if do_FitPRF_perSession
    % outputfolder
    result_folder = ['FitResult_sub-' MONKEY];
    warning off; mkdir(result_folder); warning on;
    
    % get the brain mask
    fprintf('Unpacking BrainMask');
    %uz_nii=gunzip(BrainMask_file);
    mask_nii=load_nii(BrainMask_file);%load_nii(uz_nii{1});
    %delete(uz_nii{1});
    fprintf(' ...done\n');
    
    % run the model-fits
    for s=length(sessions):-1:1
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
        end
        
        % fit pRF -----
        options.vxs = find(mask_nii.img>0);
        Sess(s).result = analyzePRF(stimulus,fmri_data,TR,options);
        
        % save the result ----
        fprintf('Saving the result: ');
        save(fullfile(result_folder,['pRF_Sess-' sessions{s}]),...
            'Sess','-v7.3');
        % also save as nifti files
        % angle ---
        fprintf('Angles ');
        nii = make_nii(Sess(s).result.ang,[1 1 1],[],[],...
            'pRF fit: Angles (deg)');
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_ang.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_ang.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_ang.nii']));
        % ecc ---
        fprintf('Ecc ');
        nii = make_nii(Sess(s).result.ecc,[1 1 1],[],[],...
            'pRF fit: Eccentricity (pix)');
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_ecc.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_ecc.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_ecc.nii']));
        % size ---
        fprintf('Size ');
        nii = make_nii(Sess(s).result.rfsize,[1 1 1],[],[],...
            'pRF fit: RF size (pix)');
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_rfsize.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_rfsize.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_rfsize.nii']));
        % R^2 Goodness of fit ---
        fprintf('R2 ');
        nii = make_nii(Sess(s).result.R2,[1 1 1],[],[],...
            'pRF fit: R2 Goodnes off fit');
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_R2.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_R2.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_R2.nii']));
        
        fprintf('>> Done!\n');
    end
end
