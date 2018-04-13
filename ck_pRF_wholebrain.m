%% WHICH DATA =============================================================
Danny_pRF_goodruns;
do_Resave = true;
do_FitPRF_perSession = true;

% This is the 'pure' sweep to volume map
SwVolMap = {    1 , 6:25    ;...
                2 , 41:60   ;...
                3 , 62:80   ;...
                4 , 96:115  ;...
                5 , 116:135 ;...
                6 , 151:170 ;...    
                7 , 171:190 ;...
                8 , 206:220 };
            
%% INITIALIZE =============================================================
% Add nifti reading toolbox
addpath(genpath('/Users/chris/Documents/MATLAB/TOOLBOX/NIfTI'))
% Add Kendrick Kay's pRF analysis toolbox
addpath(genpath('/Users/chris/Documents/MRI_ANALYSIS/analyzePRF'))
% Link to the brain mask
if strcmp(MONKEY, 'danny')
    BrainMask_file = ['/NHP-MRI/NHP-BIDS/manual-masks/sub-danny/'...
        'ses-20180117/func/T1_to_func_brainmask_zcrop.nii.gz'];
else
    error('Unknown monkey name or no mask available')
end

% create a folder to save outputs in
out_folder = ['pRF_sub-' MONKEY];
warning off
mkdir(out_folder);
warning on

%% GET THE FILE-PATHS OF THE IMAGING  & STIM-MASK FILES ===================
% All functional runs that are preprocessed with the BIDS pipeline are 
% resampled to 1x1x1 mm isotropic voxels, reoriented from sphinx,
% motion corrected, (potentially smoothed with 2 mm FWHM), and
% registered to an example functional volume 
% (so they're already in a common space)
% do the analysis in this functional space than we can register to hi-res
% anatomical data and/or the NMT template later
sessions = unique(DATA(:,1)); 
monkey_path_nii = ['/NHP-MRI/NHP-BIDS/derivatives/'...
    'featpreproc/highpassed_files/sub-' MONKEY ];
monkey_path_stim = ['/NHP-MRI/NHP-BIDS/sub-' MONKEY ];
for s=1:length(sessions)
    sess_path_nii{s} = [monkey_path_nii '/ses-' sessions{s} '/func']; %#ok<*SAGROW>
    sess_path_stim{s} = [monkey_path_stim '/ses-' sessions{s} '/func'];
    runs = unique(DATA(strcmp(DATA(:,1),sessions{s}),2));    
    for r=1:length(runs)
        run_path_nii{s,r} = [sess_path_nii{s} '/' ls([sess_path_nii{s} ...
            '/*run-' runs{r} '*'])];
        run_path_stim{s,r} = [sess_path_stim{s} '/' ls([sess_path_stim{s} ...
            '/*run-' runs{r} 'model*/StimMask.mat'])];
        sweepinc{s,r} = DATA( ...
            (strcmp(DATA(:,1),sessions{s}) & strcmp(DATA(:,2),runs{r})),3);
    end
end

%% LOAD & RE-SAVE STIMULUS MASKS & NIFTI ==================================
if do_Resave
    for s=1:size(run_path_stim,1)
        for r=1:size(run_path_stim,2)
            % stimulus mask -----
            load(run_path_stim{s,r});
            % loads variable called stimulus (x,y,t) in volumes
            sinc = sweepinc{s,r};
            firstvol = SwVolMap{min(sinc),2}(1) - 5;
            lastvol = SwVolMap{max(sinc),2}(end) + 5;
            vinc=firstvol:lastvol;
            % volumes ------
            temp_nii=load_nii(run_path_nii{s,r});
            % save the session-based stims & vols -----
            for v=1:length(vinc)
                s_run(r).stim{v} = stimulus(:,:,v);
                s_run(r).vol{v} = temp_nii(:,:,:,v);
            end
            clear stimulus temp_nii
        end
        save([out_folder '/ses-' sessions{s} ],'s_run','-v7.3');
        clear s_run
    end
else
    fprintf('Not doing the resave (assuming this has already been done)\n')
end
        
%% MODEL PRFs /SESSION ====================================================
% Do the pRF model fit on a session/day basis. Concatenate runs.
% Modeling everything together may be overkill (size wise)
if do_FitPRF_perSession
    for s=1:length(sessions)
        % load data -----
        load([out_folder '/ses-' sessions{s} ]); % s_run
        
        % concatenate -----
        stimulus={};vols = {};
        for r=1:length(s_run)
            stimulus() = s_run(r).stim
            
            
        
        % fit pRF -----
        Sess_result{s} = analyzePRF(stimulus,data,tr);
        
        % 
        
        
    end
end
