clear all; clc; %#ok<*CLALL>

Monkey = 'lick'; %#ok<*UNRCH>
MONKEY = Monkey; MONKEY(1)=upper(MONKEY(1));
Session = 'mua';
Instance = 1;
numWorkers = 2;
modeltype = 'linear_ephys';
resfld = ['ephys_xval_test'];
AllowNegGain = false;
TypicalGain = 1;
MaxIter = 1000;
SignalGain=1;

%fprintf('RUNNING IN DEBUG MODE! CHANGE THIS FLAG FOR PRODUCTION!\n');
%%
for Instance = 1:8
    fprintf(['Instance ' num2str(Instance) '\n'])
    
    %% These are fixed for this configuration =================================
    TR=0.5; % temporal resolution of the ephys signal
    mlroot = pwd; % this is $TMPDIR/PRF when running it on LISA (fast disks)
    Pix2Deg = 1/29.5;
    
    %% Prep & Load ========================================================
    
    % change characters to numbers
    if ischar(Instance); Instance = eval(Instance); end
    if ischar(numWorkers); numWorkers = eval(numWorkers); end
    InstanceLabel = num2str(Instance);
    
    % Notification of the fact that we're starting
    disp(['Starting script for job ' Monkey ', ' Session ', Instance ' ...
        num2str(Instance)])
    
    result_folder = fullfile(mlroot, 'Results', Monkey, resfld, ...
        ['Instance_' InstanceLabel]);
    
    % make outputfolder
    %result_folder = fullfile(mlroot, 'TestResults', Monkey, Session, ...
    %   ['Instance_' InstanceLabel]);
    
    if ~exist(result_folder,'dir')
        [~,~,~]=mkdir(result_folder);
    end
    %fprintf(['Saving results in: ' result_folder '\n']);
    
    %% MODEL PRFs /SESSION ================================================
    addpath(genpath('~/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/LISA/PRF/Code'));
    
    fprintf(['=== Fitting pRF model for ses-' Session ' ===\n']);
    fprintf('Loading data...\n');
    
    if ispc
        datafld = ['\\vs02\VandC\NHP_MRI\Projects\pRF\Data\ephys\' Monkey];
    else
        datafld = ['~/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/LISA/PRF/Data/ephys/' Monkey];
    end
    
    % load the stimulus
    f=dir(fullfile(datafld,'stim',[MONKEY '*']));
    load(fullfile(f.folder,f.name),'STIM');
    
    if strcmp(Session, 'mua')
        % load the responses
        if strcmp(Monkey,'lick')
            RESP0=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mMUA.mat']));
            C0=RESP.C;
        elseif strcmp(Monkey,'aston')
            RESP0=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mMUA.mat']));
            C0=RESP.C;
        end
        
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
        
        % concatenate -----
        stimulus0={};ephys_data0={};
        stimulus={};ephys_data={};
        fprintf('Concatenating stimuli and volumes...\n');
        
        ephys_data0{1} = [];
        for ch=1:128
            ephys_data0{1}=cat(1,ephys_data0{1},...
                (RESP0.mMUA(ch).bar - RESP0.mMUA(ch).BL).*SignalGain);
        end
        stimulus0{1}=[];
        for imgnr=1:length(STIM.img)
            % RESAMPLE STIMULUS >> 295 x 295 means 10px = 1 deg
            %rsIMG = imresize(STIM.img{imgnr} ,[295 295]);
            %stimulus{1}=cat(3,stimulus{1},rsIMG);
            stimulus0{1}=cat(3,stimulus0{1},STIM.img{imgnr});
        end
        
        ephys_data{1}=[];ephys_data{2}=[];
        for ch=1:128
            ephys_data{1}=cat(1,ephys_data{1},...
                (RESP.mMUA_odd(ch).bar - RESP.mMUA_odd(ch).BL)*SignalGain);
            ephys_data{2}=cat(1,ephys_data{2},...
                (RESP.mMUA_even(ch).bar - RESP.mMUA_even(ch).BL)*SignalGain);
        end
        stimulus{1}=[];stimulus{2}=[];
        for imgnr=1:length(STIM.img)
            % RESAMPLE STIMULUS >> 295 x 295 means 10px = 1 deg
            %rsIMG = imresize(STIM.img{imgnr} ,[295 295]);
            %stimulus{1}=cat(3,stimulus{1},rsIMG);
            stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
            stimulus{2}=stimulus{1};
        end
        
        %         % fit pRF -----
        %         % get indices to mask voxels > 0
        %         options.display = 'final';
        %
        %         % set crossvalidation option
        %         options.xvalmode = cv; % two-fold cross-validation (first half of runs; second half of runs)
        %
        %         % no denoising
        %         options.wantglmdenoise = 0;
        %
        %         % set typicalgain to a lower value
        %         options.typicalgain = TypicalGain;
        %
        %         % maximum number of fitting iterations
        %         options.maxiter = MaxIter;
        %
        %         % allow negative gain factors
        %         options.allowneggain = false;
        %
        %         % no drift correction for ephys
        %         options.maxpolydeg = 0;
        %
        %         % start a parallel pool of workers
        %         if ~isempty(numWorkers)
        %             parpool(numWorkers);
        %         else
        %             % if numWorkers = []
        %             % don't predefine the number of workers
        %             % let it take the max available when running
        %         end
        %
        %         % run analyzePRF tool
        %         result = analyzePRF_modeltype(stimulus,ephys_data,TR,options,modeltype);
        %         result.Chan = C;
        %         result.Pix2Deg = Pix2Deg;
        %
        %         % save the result ----
        %         fprintf('\nSaving the result: ');
        %         save(fullfile(result_folder,['pRF_Sess-' Session '_Inst_' num2str(Instance)]),...
        %             'result','-v7.3');
        %         %cd ..
        fprintf('>> Done!\n');
        
    elseif strcmp(Session,'lfp')
        % load the responses
        
        if strcmp(Monkey,'lick')
            RESP0=load(fullfile(datafld,Session,['Lick_20180807_B2_array_' num2str(Instance) '_mLFP.mat']));
            C0=RESP.C;
        elseif strcmp(Monkey,'aston')
            RESP0=load(fullfile(datafld,Session,['Aston_20181004_B1_array_' num2str(Instance) '_mLFP.mat']));
            C0=RESP.C;
        end
        
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
        
        for fb=1:5 % loop over frequency bands
            % concatenate -----
            stimulus0={};ephys_data0={};
            stimulus={};ephys_data={};
            fprintf(['Frequency band ' num2str(fb) '\n']);
            fprintf('Concatenating stimuli and volumes...\n');
            
            ephys_data0{1} = [];
            for ch=1:128
                ephys_data0{1}=cat(1,ephys_data0{1},...
                    RESP0.mLFP(ch).freq(fb).bar - ...
                    RESP0.mLFP(ch).freq(fb).BL);
            end
            stimulus0{1}=[];
            for imgnr=1:length(STIM.img)
                stimulus0{1}=cat(3,stimulus0{1},STIM.img{imgnr});
            end
            
            ephys_data{1}=[];ephys_data{2}=[];
            for ch=1:128
                ephys_data{1}=cat(1,ephys_data{1},...
                    RESP.mLFP_odd(ch).freq(fb).bar - ...
                    RESP.mLFP_odd(ch).freq(fb).BL);
                ephys_data{2}=cat(1,ephys_data{2},...
                    RESP.mLFP_even(ch).freq(fb).bar - ...
                    RESP.mLFP_even(ch).freq(fb).BL);
            end
            stimulus{1}=[];stimulus{2}=[];
            for imgnr=1:length(STIM.img)
                stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
                stimulus{2}=stimulus{1};
            end
            
            
            %             % fit pRF -----
            %             % get indices to mask voxels > 0
            %             options.display = 'final';
            %
            %             % set crossvalidation option
            %             options.xvalmode = cv; % two-fold cross-validation (first half of runs; second half of runs)
            %
            %             % no denoising
            %             options.wantglmdenoise = 0;
            %
            %             % set typicalgain to a lower value
            %             options.typicalgain = TypicalGain;
            %
            %             % maximum number of fitting iterations
            %             options.maxiter = MaxIter;
            %
            %             % allow negative gain factors
            %             options.allowneggain = AllowNegGain;
            %
            %             % no drift correction for ephys
            %             options.maxpolydeg = 0;
            %
            %             % start a parallel pool of workers
            %             if ~isempty(numWorkers)
            %                 parpool(numWorkers);
            %             else
            %                 % if numWorkers = []
            %                 % don't predefine the number of workers
            %                 % let it take the max available when running
            %             end
            %
            %             % run analyzePRF tool
            %             result = analyzePRF_modeltype(stimulus,ephys_data,TR,options,modeltype);
            %             result.Chan = C;
            %             result.Pix2Deg = Pix2Deg;
            %
            %             % save the result ----
            %             fprintf('\nSaving the result: ');
            %             save(fullfile(result_folder,['pRF_Sess-' Session '_fb' ...
            %                 num2str(fb) '_Inst_' num2str(Instance)]),'result','-v7.3');
            %             delete(gcp('nocreate'));
        end
        %cd ..
        fprintf('>> Done!\n');
    end
    delete(gcp('nocreate'));
end

%% Plot the odd and even response profiles
figure;
for i=1:1024
    subplot(32,32,i);hold on;
    plot(ephys_data{1}(1,:),'r')
    plot(ephys_data{2}(1,:),'b')
    plot(ephys_data0{1}(1,:),'k')
end

%% print some results
fprintf(['Angles: ' num2str(result.ang(1)) ' ' num2str(result.ang(2)) '\n'])
fprintf(['Ecc: ' num2str(result.ecc(1)) ' ' num2str(result.ecc(2)) '\n'])
fprintf(['R2: ' num2str(result.R2(1)) ' ' num2str(result.R2(2)) '\n'])
fprintf(['Sz: ' num2str(result.rfsize(1)) ' ' num2str(result.rfsize(2)) '\n'])
fprintf(['Gain: ' num2str(result.gain(1)) ' ' num2str(result.gain(2)) '\n'])

%% plot the prediction and data
% Perform some setup
if 1
    data = ephys_data;
    tr=TR;
    
    % Define some variables
    res = [size(stimulus{1},1) size(stimulus{1},2)];
    %res = [295 295];                    % row x column resolution of the stimuli
    resmx = max(res);                   % maximum resolution (along any dimension)
    degs = result.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model
    
    % Pre-compute cache for faster execution
    [d,xx,yy] = makegaussian2d(resmx,2,2,2,2);
    
    % Prepare the stimuli for use in the model
    stimulusPP = {};
    for p=1:length(stimulus)
        stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
        stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
    end
    
    % Define the model function.  This function takes parameters and stimuli as input and
    % returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
    % of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
    % Although it looks complex, what the function does is pretty straightforward: construct a
    % 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
    % Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
    % taking care to not bleed over run boundaries.
    %     modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
    %         makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
    %         (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5));
    modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
        makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
        (2*pi*abs(pp(3))^2))); 0]);
    
    % Construct projection matrices that fit and remove the polynomials.
    % Note that a separate projection matrix is constructed for each run.
    polymatrix = {};
    for p=1:length(degs)
        polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
    end
    
    % Inspect the data and the model fit
    
    %% Which channel should we inspect?
    ch = 1;
    
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {};
    for p=1:length(data)
        datats{p} =  polymatrix{p}*data{p}(ch,:)';
        modelts{p} = polymatrix{p}*modelfun(result.params(1,:,ch),single(stimulusPP{p}));
    end
    
    % Visualize the results
    figure; hold on;
    %set(gcf,'Units','points','Position',[100 100 1000 100]);
    plot(cat(1,datats{:}),'r-');
    plot(cat(1,modelts{:}),'b-');
    straightline(240,'v','k-');
    xlabel('Time (s)');
    ylabel('Signal');
    ax = axis;
    %axis([.5 1200+.5 ax(3:4)]);
    title('Time-series data');
    
    
end




