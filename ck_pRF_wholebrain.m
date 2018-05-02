function ck_pRF_wholebrain(SessionList, do_Resave, do_FitPRF_perSession)
% collects data, concatenates, downsamples stimulus and resaves
% fits the prf model to voxels

doUpsample = true; % upsamples to twice TR, so 1.25s
doExtraRegression = true; % include motion information as regressor
fitOnlyPosterior = true; % mask out anterior part of the brain to speed things up
numWorkers = [];    % set the number of parallel processes
                    % [] is default max number available
                    % sometimes you'd want less because of memory limitations

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

%% Volume map
TR=2.5;

% This is the 'pure' sweep to volume map
% For the 230 volume version (15 vol blanks)
SwVolMap_new = { ...
    1 , 6:25    ;...
    2 , 41:60   ;...
    3 , 61:80   ;...
    4 , 96:115  ;...
    5 , 116:135 ;...
    6 , 151:170 ;...
    7 , 171:190 ;...
    8 , 206:225 };
% For the 210 volume version (10 vol blanks)
SwVolMap_old = { ...
    1 , 6:25    ;...
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
if doUpsample
    out_folder = ['pRF_sub-' MONKEY '_us'];
else
    out_folder = ['pRF_sub-' MONKEY]; %#ok<*UNRCH>
end
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
sessions = unique(DATA(:,1)); %#ok<*NODEF>

monkey_path_nii = fullfile(BIDS_basepath, 'derivatives',...
    'featpreproc','highpassed_files',['sub-' MONKEY]);
monkey_path_stim = fullfile(BIDS_basepath,['sub-' MONKEY]);

monkey_path_motion.regress = fullfile(BIDS_basepath, 'derivatives',...
    'featpreproc','motion_corrected',['sub-' MONKEY]);
monkey_path_motion.outlier = fullfile(BIDS_basepath, 'derivatives',...
    'featpreproc','motion_outliers',['sub-' MONKEY]);

for s=1:length(sessions)
    sess_path_nii{s} = fullfile(monkey_path_nii, ['ses-' sessions{s}], 'func'); %#ok<*SAGROW>
    sess_path_stim{s} = fullfile(monkey_path_stim, ['ses-' sessions{s}], 'func');
    sess_path_motreg{s} = fullfile(monkey_path_motion.regress, ['ses-' sessions{s}], 'func');
    sess_path_motout{s} = fullfile(monkey_path_motion.outlier, ['ses-' sessions{s}], 'func');
    runs = unique(DATA(strcmp(DATA(:,1),sessions{s}),2));
    for r=1:length(runs)
        if ispc % the ls command works differently in windows
            a=ls( fullfile(sess_path_nii{s},['*run-' runs{r} '*.nii.gz']));
            run_path_nii{s,r} = fullfile(sess_path_nii{s},a(1:end-3));
            b = ls( fullfile(sess_path_stim{s}, ...
                ['*run-' runs{r} '*model*']));
            run_path_stim{s,r}= fullfile(sess_path_stim{s},b,'StimMask.mat');
            c=ls( fullfile(sess_path_motreg{s},['*run-' runs{r} '*.param.1D']));
            run_path_motreg{s,r} = fullfile(sess_path_motreg{s},c);
            d=ls( fullfile(sess_path_motout{s},['*run-' runs{r} '*.outliers.txt']));
            run_path_motout{s,r} = fullfile(sess_path_motout{s},d);
            e = ls( fullfile(sess_path_stim{s}, ...
                ['*run-' runs{r} '*model*']));
            run_path_rew{s,r}= fullfile(sess_path_stim{s},e,'RewardEvents.txt');
        else
            a = ls( fullfile(sess_path_nii{s},['*run-' runs{r} '*.nii.gz']));
            run_path_nii{s,r} = a(1:end-3);
            run_path_nii{s,r} = run_path_nii{s,r}(1:end-1);
            run_path_stim{s,r} = ls( fullfile(sess_path_stim{s}, ...
                ['*run-' runs{r} '*model*'],'StimMask.mat'));
            run_path_stim{s,r} = run_path_stim{s,r}(1:end-1);
            run_path_motreg{s,r} = ls( fullfile(sess_path_motreg{s}, ...
                ['*run-' runs{r} '*.param.1D']));
            run_path_motreg{s,r}=run_path_motreg{s,r}(1:end-1);
            run_path_motout{s,r} = ls( fullfile(sess_path_motout{s}, ...
                ['*run-' runs{r} '*_outliers.txt']));
            run_path_motout{s,r}=run_path_motout{s,r}(1:end-1);
            run_path_rew{s,r} = ls( fullfile(sess_path_stim{s}, ...
                ['*run-' runs{r} '*model*'],'RewardEvents.txt'));
            run_path_rew{s,r}=run_path_rew{s,r}(1:end-1);
        end
        sweepinc{s,r} = DATA( ...
            (strcmp(DATA(:,1),sessions{s}) & strcmp(DATA(:,2),runs{r})),3);
    end
end

%% LOAD & RE-SAVE STIMULUS MASKS & NIFTI ==================================
if do_Resave
    for s=1:size(run_path_stim,1) % sessions
        fprintf(['Processing session ' sessions{s} '\n']);
        rps = [];
        for i=1:size(run_path_stim,2) 
            if ~isempty(run_path_stim{s,i}) 
                rps=[rps i]; 
            end 
        end
        for r=rps % runs
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
            
            % motion regressors
            if doExtraRegression
                fprintf('Processing motion regressors\n');
                
                % outliers ---
                sss=dir(run_path_motout{s,r}); fsz=sss.bytes;
                if fsz>0 % file not empty
                    s_run(r).motion.outliers = dlmread(run_path_motout{s,r});
                    % these volumes can potentially be removed 
                else
                    s_run(r).motion.outliers = []; % no outliers
                end

                % translation/rotation ---
                % 1) x-shift  2) y-shift  
                % 3) z-shift  4) z-angle  
                % 5) x-angle  6) y-angle  
                % 7) x-scale  8) y-scale  
                % 9) z-scale  10) y/x-shear  
                % 11) z/x-shear  12) z/y-shear 
                motion.estimates=dlmread(run_path_motreg{s,r},'',2,0);
                motion.estimates=motion.estimates(vinc,:);
                s_run(r).motion.estimates=motion.estimates-motion.estimates(1,:);
                
                % reward events ---
                sss=dir(run_path_rew{s,r}); fsz=sss.bytes;
                if fsz>0 % file not empty
                    rew_ev = dlmread(run_path_rew{s,r});
                else
                    rew_ev = []; % no outliers
                end
                % convert to reward per volume (based on TR)
                rew_reg =[];
                rew_ev(:,4) = rew_ev(:,1)+rew_ev(:,2); % reward end moment
                for vv=vinc
                    tw = [(vv*TR)-TR vv*TR];
                    ind=rew_ev(:,1)>tw(1) & rew_ev(:,1)<tw(2);
                    if sum(ind)>0
                        rew_reg = [rew_reg; sum(rew_ev(ind,2))];
                    else
                        rew_reg = [rew_reg; 0];
                    end
                end
                s_run(r).rew = rew_reg;
            end
            
            % save the session-based stims & vols -----
            for v=1:length(vinc)
                % resample image (160x160 pix gives 10 pix/deg)
                s_run(r).stim{v} = imresize(stimulus(:,:,vinc(v)),[160 160]);
                s_run(r).vol{v} = temp_nii.img(:,:,:,vinc(v));
            end
            clear stimulus temp_nii

            
            % if requested, upsample temporal resolution
            if doUpsample
                % stim ---
                tempstim = s_run(r).stim;
                ups_stim = cell(1,2*length(tempstim));
                ups_stim(1:2:end) = tempstim; 
                ups_stim(2:2:end) = tempstim;
                s_run(r).stim = ups_stim;
                clear tempstim ups_stim
                
                % bold ---
                us_nii=[];
                for v=1:length(s_run(r).vol)
                    us_nii=cat(4,us_nii,s_run(r).vol{v});
                end
                fprintf('Upsampling BOLD data...\n');
                us_nii = tseriesinterp(us_nii,TR,TR/2,4);
                for v=1:size(us_nii,4)
                    s_run(r).vol{v} = us_nii(:,:,:,v);
                end
                clear us_nii
                
                % motion regressors ---
                if doExtraRegression
                    tempmot = s_run(r).motion.estimates;
                    ups_mot = nan(size(tempmot,1)*2,size(tempmot,2));
                    ups_mot(1:2:end,:) = tempmot; 
                    ups_mot(2:2:end,:) = tempmot;
                    s_run(r).motion.estimates = ups_mot;
                    clear tempstim ups_mot
                    
                    temprew = s_run(r).rew;
                    ups_rew = nan(size(temprew,1)*2,size(temprew,2));
                    ups_rew(1:2:end,:) = temprew; 
                    ups_rew(2:2:end,:) = temprew;
                    s_run(r).rew = ups_rew;
                    clear tempstim ups_mot
                end
            end
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
    
    % get the brain mask
    fprintf('Unpacking BrainMask');
    %uz_nii=gunzip(BrainMask_file);
    mask_nii=load_nii(BrainMask_file);%load_nii(uz_nii{1});
    if fitOnlyPosterior
        mask_nii.img(:,50:end,:)=0; % mask out anterior part (speed up fits)
    end
    
    %delete(uz_nii{1});
    fprintf(' ...done\n');
    
    % run the model-fits
    
    % 1 - ses-20171116.mat
    % 2 - ses-20171129.mat
    % 3 - ses-20171207.mat
    % 4 - ses-20171214.mat
    % 5 - ses-20171220.mat
    % 6 - ses-20180117.mat
    % 7 - ses-20180124.mat
    % 8 - ses-20180125.mat
    % 9 - ses-20180131.mat
    % 10 - ses-20180201.mat
    
    session_order = 1;%[2 4:10]; % you can do the bad ones later (or not)
    for s=session_order %length(sessions):-1:1
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
        if ~isempty(numWorkers)
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
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_ang.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_ang.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_ang.nii']));
        % ecc ---
        fprintf('Ecc ');
        nii = make_nii(result.ecc,[1 1 1],[],[],...
            'pRF fit: Eccentricity (pix)');
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_ecc.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_ecc.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_ecc.nii']));
        % size ---
        fprintf('Size ');
        nii = make_nii(result.rfsize,[1 1 1],[],[],...
            'pRF fit: RF size (pix)');
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_rfsize.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_rfsize.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_rfsize.nii']));
        % R^2 Goodness of fit ---
        fprintf('R2 ');
        nii = make_nii(result.R2,[1 1 1],[],[],...
            'pRF fit: R2 Goodnes off fit');
        save_nii(nii, fullfile(result_folder, ['Sess-' sessions{s}' '_R2.nii']));
        gzip(fullfile(result_folder, ['Sess-' sessions{s}' '_R2.nii']));
        delete(fullfile(result_folder, ['Sess-' sessions{s}' '_R2.nii']));
        
        fprintf('>> Done!\n');
    end
end
