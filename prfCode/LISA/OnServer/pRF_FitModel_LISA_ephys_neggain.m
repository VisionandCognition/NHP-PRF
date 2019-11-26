function pRF_FitModel_LISA_ephys_neggain(Monkey,Session,Instance,numWorkers,modeltype,cv,resfld)
% This is the script that runs the pRF fit
% it should be compiled and called from run_CompiledMatlab_LISA.sh
% This means it cannot use 'addpath', instead toolboxes or required
% function should be added when it is compiled usin -a /toolbox etc.

% compile on LISA with:
% module load matlab
% echo "mcc -m $HOME/PRF/Code/pRF_FitModel_LISA.m -a $HOME/PRF/Code/analyzePRF " | matlab -nodisplay
% module unload matlab

% The uncompiled version of this function should be in
% $HOME/PRF/Code

% this compiled function will be running from:
% "$TMPDIR"/PRF

% Data will be in
% "$TMPDIR"/PRF

DEBUG=false;
if DEBUG
    Monkey = 'lick'; %#ok<*UNRCH>
    Session = 'mua';
    Instance = 1;
    numWorkers = 2;
    modeltype = 'linear_ephys';
    cv = 1;
    resfld = 'debug'
    fprintf('RUNNING IN DEBUG MODE! CHANGE THIS FLAG FOR PRODUCTION!\n');
end

%% These are fixed for this configuration =================================
TR=0.5;
mlroot = pwd; % this is $TMPDIR/PRF when running it on LISA (fast disks)
Pix2Deg = 1/29.5;

%% Prep & Load ============================================================

% change characters to numbers
if ischar(Instance); Instance = eval(Instance); end
if ischar(numWorkers); numWorkers = eval(numWorkers); end
if ischar(cv); cv = eval(cv); end

InstanceLabel = num2str(Instance);

% Notification of the fact that we're starting
disp(['Starting script for job ' Monkey ', ' Session ', Instance ' num2str(Instance)])

% make outputfolder
% result_folder = [mlroot '/Results/' Monkey '/' Session '/Slices_' SliceLabel];
result_folder = fullfile(mlroot, 'Results', Monkey, resfld, ['Instance_' InstanceLabel]);

if ~exist(result_folder,'dir')
    [~,~] = mkdir(result_folder);
end
fprintf(['Saving results in: ' result_folder '\n']);

%% MODEL PRFs /SESSION ====================================================
fprintf(['=== Fitting pRF model for ses-' Session ' ===\n']);
fprintf('Loading data...\n');

datafld = fullfile('/home/pcklink/PRF/Data/ephys/', Monkey);

% load the stimulus
if strcmp(Monkey,'lick')
    load(fullfile(datafld,'stim','Lick_20180807_B2_STIM.mat'))
elseif strcmp(Monkey,'aston')
    load(fullfile(datafld,'stim','Aston_20181004_B1_STIM.mat'))
end

if strcmp(Session, 'mua')
    % load the responses
    if cv == 0 % no crossvalidation
        if strcmp(Monkey,'lick')
            RESP=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mMUA.mat']));
            C=RESP.C;
        elseif strcmp(Monkey,'aston')
            RESP=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mMUA.mat']));
            C=RESP.C;
        end
    elseif cv == 1 % crossvalidation
        if strcmp(Monkey,'lick')
            RESP=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mMUA_odd.mat']));
            C=RESP.C;
            RESP2=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mMUA_even.mat']));
            RESP.mMUA_even=RESP2.mMUA_even;
        elseif strcmp(Monkey,'aston')
            RESP=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mMUA_odd.mat']));
            C=RESP.C;
            RESP2=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mMUA_even.mat']));
            RESP.mMUA_even=RESP2.mMUA_even;
        end
    end
    
    % concatenate -----
    stimulus={};ephys_data={};
    fprintf('Concatenating stimuli and volumes...\n');
    
    if cv == 0
        ephys_data{1} = [];
        for ch=1:128
            %ephys_data{1}=cat(1,ephys_data{1},...
            %    RESP.mMUA(ch).bar - RESP.mMUA(ch).BL);
            ephys_data{1}=cat(1,ephys_data{1},...
                    RESP.mMUA(ch).bar.*SignalGain);
        end
        stimulus{1}=[];
        for imgnr=1:length(STIM.img)
            stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
        end
    elseif cv == 1
        ephys_data{1}=[];ephys_data{2}=[];
        for ch=1:128
%             ephys_data{1}=cat(1,ephys_data{1},...
%                 RESP.mMUA_odd(ch).bar - RESP.mMUA_odd(ch).BL);
%             ephys_data{2}=cat(1,ephys_data{2},...
%                 RESP.mMUA_even(ch).bar - RESP.mMUA_even(ch).BL);
            ephys_data{1}=cat(1,ephys_data{1},...
                    RESP.mMUA_odd(ch).bar.*SignalGain);
            ephys_data{2}=cat(1,ephys_data{2},...
                    RESP.mMUA_even(ch).bar.*SignalGain);
        end
        stimulus{1}=[];stimulus{2}=[];
        for imgnr=1:length(STIM.img)/2
            stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
            stimulus{2}=stimulus{1};
        end
        
        inv_idx = [150:-1:121 180:-1:151 210:-1:181 240:-1:211];
        ephys_data2{1}=[]; ephys_data2{2}=[];
        for elec = 1:size(ephys_data{1},1)
            ephys_data2{1} = cat(1,ephys_data2{1},...
                mean([ephys_data{1}(elec,1:120);ephys_data{1}(elec,inv_idx)],1));
            ephys_data2{2} = cat(1,ephys_data2{2},...
                mean([ephys_data{2}(elec,1:120);ephys_data{2}(elec,inv_idx)],1));
        end
        ephys_data_org = ephys_data;
        ephys_data = ephys_data2; clear ephys_data2
        
    end
    
    % fit pRF -----
    options.display = 'final';
    
    % set crossvalidation option
    options.xvalmode = cv; % two-fold cross-validation (first half of runs; second half of runs)
    
    % no denoising
    options.wantglmdenoise = 0;
    
    % set typicalgain to a lower value
    %options.typicalgain = 10;
    
    % allow negative gain factors
    options.allowneggain = true;

    % no drift correction for ephys
    options.maxpolydeg = 0;
    
    % start a parallel pool of workers
    if ~isempty(numWorkers)
        parpool(numWorkers);
    else
        % if numWorkers = []
        % don't predefine the number of workers
        % let it take the max available when running
    end
    
    % run analyzePRF tool
    result = analyzePRF_modeltype(stimulus,ephys_data,TR,options,modeltype);
    result.Chan = C;
    result.Pix2Deg = Pix2Deg;
    
    % save the result ----
    fprintf('\nSaving the result: ');
    save(fullfile(result_folder,['pRF_Sess-' Session '_Inst_' num2str(Instance)]),...
        'result','-v7.3');
    cd ..
    fprintf('>> Done!\n');
    
elseif strcmp(Session,'lfp')
    % load the responses
    if cv == 0 % no crossvalidation
        if strcmp(Monkey,'lick')
            RESP=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mLFP.mat']));
            C=RESP.C;
        elseif strcmp(Monkey,'aston')
            RESP=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mLFP.mat']));
            C=RESP.C;
        end        
    elseif cv == 1 % crossvalidation
        if strcmp(Monkey,'lick')
            RESP=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mLFP_odd.mat']));
            C=RESP.C;
            RESP2=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mLFP_even.mat']));
            RESP.mLFP_even=RESP2.mLFP_even;
        elseif strcmp(Monkey,'aston')
            RESP=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mLFP_odd.mat']));
            C=RESP.C;
            RESP2=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mLFP_even.mat']));
            RESP.mLFP_even=RESP2.mLFP_even;
        end        
    end
    
    for fb=1:5 % loop over frequency bands
        % concatenate -----
        stimulus={};ephys_data={};
        fprintf(['Frequency band ' num2str(fb) '\n']);
        fprintf('Concatenating stimuli and volumes...\n');
        
        if cv == 0
            ephys_data{1} = [];
            for ch=1:128
                ephys_data{1}=cat(1,ephys_data{1}, ...
                    RESP.mLFP(ch).freq(fb).bar - ...
                    RESP.mLFP(ch).freq(fb).BL);
            end
            stimulus{1}=[];
            for imgnr=1:length(STIM.img)
                stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
            end
        elseif cv == 1
            ephys_data{1}=[];ephys_data{2}=[];
            for ch=1:128
                ephys_data{1}=cat(1,ephys_data{1}, ...
                    RESP.mLFP_odd(ch).freq(fb).bar - ...
                    RESP.mLFP_odd(ch).freq(fb).BL);
                ephys_data{2}=cat(1,ephys_data{2},...
                    RESP.mLFP_even(ch).freq(fb).bar - ...
                    RESP.mLFP_even(ch).freq(fb).BL);
            end
            stimulus{1}=[];stimulus{2}=[];
            for imgnr=1:length(STIM.img)/2
                stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
                stimulus{2}=stimulus{1};
            end
            
            inv_idx = [150:-1:121 180:-1:151 210:-1:181 240:-1:211];
            ephys_data2{1}=[]; ephys_data2{2}=[];
            for elec = 1:size(ephys_data{1},1)
                ephys_data2{1} = cat(1,ephys_data2{1},...
                    mean([ephys_data{1}(elec,1:120);ephys_data{1}(elec,inv_idx)],1));
                ephys_data2{2} = cat(1,ephys_data2{2},...
                    mean([ephys_data{2}(elec,1:120);ephys_data{2}(elec,inv_idx)],1));
            end
            
            ephys_data_org = ephys_data;
            ephys_data = ephys_data2; clear ephys_data2
            
        end
        
        % fit pRF -----
        options.display = 'final';
        
        % set crossvalidation option
        options.xvalmode = cv; % two-fold cross-validation (first half of runs; second half of runs)
        
        % no denoising
        options.wantglmdenoise = 0;
        
        % set typicalgain to a lower value
        %options.typicalgain = 10;
        
        % allow negative gain factors
        options.allowneggain = true;
    
        % no drift correction for ephys
        options.maxpolydeg = 0;
        
        % start a parallel pool of workers
        if ~isempty(numWorkers)
            parpool(numWorkers);
        else
            % if numWorkers = []
            % don't predefine the number of workers
            % let it take the max available when running
        end
        
        % run analyzePRF tool
        result = analyzePRF_modeltype(stimulus,ephys_data,TR,options,modeltype);
        result.Chan = C;
        result.Pix2Deg = Pix2Deg;
        
        % save the result ----
        fprintf('\nSaving the result: ');
        save(fullfile(result_folder,['pRF_Sess-' Session '_fb' ...
            num2str(fb) '_Inst_' num2str(Instance)]),'result','-v7.3');
    end
    cd ..
    fprintf('>> Done!\n');
end