function pRF_FitModel_LISA_cv(Monkey,Session,Slices,HRF,numWorkers,modeltype,cv,resfld)
% This is the script that runs the pRF fit
% it should be compiled and called from run_CompiledMatlab_LISA.sh
% This means it cannot use 'addpath', instead toolboxes or required
% function should be added when it is compiled usin -a /toolbox etc.

% compile on LISA with:
% module load matlab
% echo "mcc -m $HOME/PRF/Code/pRF_FitModel_LISA.m -a $HOME/PRF/Code/HRF -a $HOME/PRF/Code/analyzePRF -a $HOME/PRF/Code/NIfTI" | matlab -nodisplay
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
% <Monkey> string, no caps
% <Session> string YYYYMMDD
% <Slices> range of slices x:xx
% <HRF> HRF function to use (won't be used if <modeltype> is an ephys variant
% <numWorkers> can be used to cap the # of parallel processes, leave [] to
% maximize
% <modeltype> specifies what kind of model to fit ('css','css_ephys','linear','linear_ephys')
% <cv> type of crossvalidation (0,1,2)

DEBUG=false;
if DEBUG
    Monkey = 'eddy'; %#ok<*UNRCH>
    Session = 'AllSessions-avg-cv';
    Slices = 01:5;
    HRF = 'HRF_monkey';
    numWorkers = [];
    modeltype = 'css';
    cv = 0;
    fprintf('RUNNING IN DEBUG MODE! CHANGE THIS FLAG FOR PRODUCTION!\n');
end

%% These are fixed for this configuration =================================
TR=2.5;
if strcmp(Monkey, 'eddy') && strcmp(Session,'20160721')
    TR=3;
end

doUpsample=true;
mlroot = pwd; % this is $TMPDIR/PRF when running it on LISA (fast disks)

%% Prep & Load ============================================================

% change characters to numbers
if ischar(Slices); slices = eval(Slices); else slices = Slices; end
if ischar(numWorkers); numWorkers = eval(numWorkers); end
if ischar(cv); cv = eval(cv); end

SliceLabel = [num2str(slices(1),'%02d') '-' num2str(slices(end),'%02d')];

% Notification of the fact that we're starting
disp(['Starting script for job ' Monkey ', ' Session ', Slices ' Slices])

% Link to the brain mask
if strcmp(Monkey, 'danny')
    BrainMask_file = [mlroot '/T1_to_func_brainmask_zcrop.nii'];
elseif strcmp(Monkey, 'eddy')
    BrainMask_file = [mlroot '/HiRes_to_T1_mean.nii_shadowreg_Eddy_brainmask.nii'];
end

% make outputfolder
% result_folder = [mlroot '/Results/' Monkey '/' Session '/Slices_' SliceLabel];
result_folder = [mlroot '/Results/' Monkey '/' resfld '/Slices_' SliceLabel];

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

fprintf(['=== Fitting pRF model for ses-' Session ' ===\n']);
fprintf('Loading data...\n');

load(fullfile(mlroot,Session)); %#ok<*LOAD>
% The crucial variables in thes files are:
% stim  >> a structure with fields norm and inv
%       norm >> stimulus in normal order. cell array. 
%       inv >> stimulus in inverse order. cell array. 
% sess_wmeanBOLD >> weighted average BOLD response
%                >> 4d array of epi-volumes in time (4th dimension) 
% sess_wmeanBOLD >> weighted average BOLD response for inverse stimuli
%                >> 4d array of epi-volumes in time (4th dimension) 

% concatenate -----
stimulus={};fmri_data={};
fprintf('Concatenating stimuli and volumes...\n');

if iscell(sess_wmeanBOLD)
    for RUNNR = 1:length(sess_wmeanBOLD)
        fmri_data{(RUNNR-1)+RUNNR}=sess_wmeanBOLD{RUNNR}(:,:,slices,:);  %#ok<*IDISVAR,*NODEF>
        stimulus{(RUNNR-1)+RUNNR}=[]; 
        for voln = 1:size(stim.norm{RUNNR},2)
            stimulus{(RUNNR-1)+RUNNR} = cat(3, stimulus{(RUNNR-1)+RUNNR}, stim.norm{RUNNR}{voln}); %#ok<*AGROW>
        end
        if exist('sess_wmeanBOLD_inv') && ~isempty(sess_wmeanBOLD_inv{RUNNR}) %#ok<*EXIST>
            fmri_data{RUNNR+RUNNR}=sess_wmeanBOLD_inv{RUNNR}(:,:,slices,:);
            stimulus{RUNNR+RUNNR}=[];
            for voln = 1:size(stim.inv{RUNNR},2)
                stimulus{RUNNR+RUNNR} = cat(3, stimulus{RUNNR+RUNNR}, stim.inv{RUNNR}{voln});
            end      
        end
    end
else % for the xval = 0 data 
    fmri_data{1}=sess_wmeanBOLD(:,:,slices,:);
    stimulus{1}=[]; 
    for voln = 1:size(stim.norm,2)
        stimulus{1} = cat(3, stimulus{1}, stim.norm{voln}); %#ok<*AGROW>
    end

    if exist('sess_wmeanBOLD_inv') && ~isempty(sess_wmeanBOLD_inv) %#ok<*EXIST>
        fmri_data{2}=sess_wmeanBOLD_inv(:,:,slices,:);
        stimulus{2}=[];
        for voln = 1:size(stim.inv,2)
            stimulus{2} = cat(3, stimulus{2}, stim.inv{voln});
        end
    end    
end


% clear empty cell for non-existing inverse stimuli
if length(fmri_data)>1 && isempty(fmri_data{2})
    fmri_data(2)=[];
    stimulus(2)=[];
end

% fit pRF -----
% get indices to mask voxels > 0
options.vxs = find(mask_nii.img(:,:,slices)>0);
options.display = 'final';

if strcmp(HRF,'HRF_monkey')
    load(fullfile(mlroot,'HRF',HRF),'customHRF_rs','customHRF_rsus');
    if doUpsample 
        options.hrf = customHRF_rsus;
    else
        options.hrf = customHRF_rs;
    end
end

% set crossvalidation option
options.xvalmode = cv; % two-fold cross-validation (first half of runs; second half of runs)

% no denoising for average BOLD traces
options.wantglmdenoise = 0;

% set typicalgain to a lower value because it's %BOLD change now and not raw data
options.typicalgain = 0.25;

% allow negative gain factors
options.allowneggain = true;

% start a parallel pool of workers
if ~isempty(numWorkers)
    parpool(numWorkers);
else
    % if numWorkers = []
    % don't predefine the number of workers
    % let it take the max available when running
end

% run analyzePRF tool
if doUpsample % tr = TR/2
    result = analyzePRF_modeltype(stimulus,fmri_data,TR/2,options,modeltype);
else
    result = analyzePRF_modeltype(stimulus,fmri_data,TR,options,modeltype);
end

load(fullfile(mlroot,[Monkey '_refhdr']));
result.hdr_ref = ref_hdr;

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
save_nii(nii, fullfile(result_folder, [Session '_ang.nii']));
gzip(fullfile(result_folder, [Session '_ang.nii']));
delete(fullfile(result_folder, [Session '_ang.nii']));
% ecc ---
fprintf('Ecc ');
nii = make_nii(result.ecc,[1 1 1],[],32,...
    'pRF fit: Eccentricity (pix)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(result_folder, [Session '_ecc.nii']));
gzip(fullfile(result_folder, [Session '_ecc.nii']));
delete(fullfile(result_folder, [Session '_ecc.nii']));
% size ---
fprintf('Size ');
nii = make_nii(result.rfsize,[1 1 1],[],32,...
    'pRF fit: RF size (pix)');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(result_folder, [Session '_rfsize.nii']));
gzip(fullfile(result_folder, [Session '_rfsize.nii']));
delete(fullfile(result_folder, [Session '_rfsize.nii']));
% R^2 Goodness of fit ---
fprintf('R2 ');
nii = make_nii(result.R2,[1 1 1],[],32,...
    'pRF fit: R2 Goodnes off fit');
nii.hdr.hist = result.hdr_ref.hist;
save_nii(nii, fullfile(result_folder, [Session '_R2.nii']));
gzip(fullfile(result_folder, [Session '_R2.nii']));
delete(fullfile(result_folder, [Session '_R2.nii']));

cd ..
fprintf('>> Done!\n');