function pRF_prepdata(SessionList,doUpsample,doExtraRegression)
% collects data, concatenates, downsamples stimulus and resaves
% NB: won't run on macbook unless there's a USB-drive with all data 
%% WHICH DATA =============================================================
%clear all; clc;
if nargin <2
    fprintf('Not enough arguments specified, will use defaults:\n');
    fprintf('SessionList: Danny_pRF_goodruns\n');
    fprintf('doUpsample = true; doExtraRegression = true\n');
    pRF_PrepDatalist_Danny;
    doUpsample = true;
    doExtraRegression = true;
else
    eval(SessionList);
end

%% Sweep to volume mapping ------------------------------------------------
TR = 2.5;
TR_3 = 3.0; % 20160721

SwVolMap_230 = { ...
    1 , 6:25    ;...
    2 , 41:60   ;...
    3 , 61:80   ;...
    4 , 96:115  ;...
    5 , 116:135 ;...
    6 , 151:170 ;...
    7 , 171:190 ;...
    8 , 206:225 };

SwVolMap_218 = { ...
    1 , 6:25    ;...
    2 , 38:57   ;...
    3 , 58:77   ;...
    4 , 90:109  ;...
    5 , 110:129 ;...
    6 , 142:161 ;...
    7 , 162:182 ;...
    8 , 194:213 };

SwVolMap_210 = { ...        % also for 215 vols
    1 , 6:25    ;...
    2 , 36:55   ;...
    3 , 56:75   ;...
    4 , 86:105  ;...
    5 , 106:125 ;...
    6 , 136:155 ;...
    7 , 156:175 ;...
    8 , 186:205 };

SwVolMap_436 = { ...
    1 , 6:25    ;...
    2 , 38:57   ;...
    3 , 58:77   ;...
    4 , 90:109  ;...
    5 , 110:129 ;...
    6 , 142:161 ;...
    7 , 162:182 ;...
    8 , 194:213 ;...
    9 , 224:243 ;...
    2 , 246:265 ;...
    3 , 276:285 ;...
    4 , 308:327 ;...
    5 , 328:347 ;...
    6 , 360:379 ;...
    7 , 380:399 ;...
    8 , 412:431 };


%% INITIALIZE =============================================================
% Platform specific basepath
if ispc
    tool_basepath = 'D:\CK\code\MATLAB';
    BIDS_basepath = '\\vcnin\NHP_MRI\NHP-BIDS';
else
    tool_basepath = '~/Dropbox/MATLAB_NONGIT/TOOLBOX';
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
elseif strcmp(MONKEY, 'eddy')
    BrainMask_file = fullfile(BIDS_basepath, 'manual-masks','sub-eddy',...
        'ses-20170607b','anat','HiRes_to_T1_mean.nii_shadowreg_Eddy_brainmask.nii');
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
    sess_path_nii{s} = fullfile(monkey_path_nii, ['ses-' sessions{s}(1:8)], 'func'); %#ok<*SAGROW>
    sess_path_stim{s} = fullfile(monkey_path_stim, ['ses-' sessions{s}(1:8)], 'func');
    sess_path_motreg{s} = fullfile(monkey_path_motion.regress, ['ses-' sessions{s}(1:8)], 'func');
    sess_path_motout{s} = fullfile(monkey_path_motion.outlier, ['ses-' sessions{s}(1:8)], 'func');
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
for s=1:size(run_path_stim,1) % sessions
    fprintf(['Processing session ' sessions{s} '\n']);
    rps = [];
    if strcmp(sessions{s},'20160721')
        TR=TR_3;
    end

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
        if size(stimulus,3) == 210 || size(stimulus,3) == 215
            SwVolMap = SwVolMap_210;
        elseif size(stimulus,3) == 218
            SwVolMap = SwVolMap_218;
        elseif size(stimulus,3) == 230
            SwVolMap = SwVolMap_230;
        elseif size(stimulus,3) == 436
            SwVolMap = SwVolMap_436;           
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
            if v==1
                s_run(r).hdr = temp_nii.hdr;
            end
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