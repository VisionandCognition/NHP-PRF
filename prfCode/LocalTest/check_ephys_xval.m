clear all; clc; %#ok<*CLALL>

Monkey = 'lick'; %#ok<*UNRCH>
MONKEY = Monkey; MONKEY(1)=upper(MONKEY(1));
Session = 'mua';
resfld = ['ephys_xval_test'];
AllowNegGain = false;
TypicalGain = 1;
MaxIter = 1000;
SignalGain=1;
numWorkers=2;

%%
ephys_data0{1}=[];
ephys_data{1}=[];ephys_data{2}=[];

for Instance = 1%:8
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
            C0=RESP0.C;
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
        stimulus0={};%ephys_data0={};
        stimulus={};%ephys_data={};
        fprintf('Concatenating stimuli and data...\n');
        
        %ephys_data0{1} = [];
        for ch=1:128
            %ephys_data0{1}=cat(1,ephys_data0{1},...
            %    (RESP0.mMUA(ch).bar - RESP0.mMUA(ch).BL).*SignalGain);
            ephys_data0{1}=cat(1,ephys_data0{1},...
                (RESP0.mMUA(ch).bar).*SignalGain);
        end
        stimulus0{1}=[];
        for imgnr=1:length(STIM.img)/2
            % RESAMPLE STIMULUS >> 295 x 295 means 10px = 1 deg
            %rsIMG = imresize(STIM.img{imgnr} ,[295 295]);
            %stimulus{1}=cat(3,stimulus{1},rsIMG);
            stimulus0{1}=cat(3,stimulus0{1},STIM.img{imgnr});
        end
        
        %ephys_data{1}=[];ephys_data{2}=[];
        for ch=1:128
%             ephys_data{1}=cat(1,ephys_data{1},...
%                 (RESP.mMUA_odd(ch).bar - RESP.mMUA_odd(ch).BL)*SignalGain);
%             ephys_data{2}=cat(1,ephys_data{2},...
%                 (RESP.mMUA_even(ch).bar - RESP.mMUA_even(ch).BL)*SignalGain);
            ephys_data{1}=cat(1,ephys_data{1},...
                (RESP.mMUA_odd(ch).bar)*SignalGain);
            ephys_data{2}=cat(1,ephys_data{2},...
                (RESP.mMUA_even(ch).bar)*SignalGain);
        end
        stimulus{1}=[];stimulus{2}=[];
        for imgnr=1:length(STIM.img)
            % RESAMPLE STIMULUS >> 295 x 295 means 10px = 1 deg
            %rsIMG = imresize(STIM.img{imgnr} ,[295 295]);
            %stimulus{1}=cat(3,stimulus{1},rsIMG);
            stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
            stimulus{2}=stimulus{1};
        end

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

        end
        %cd ..
        fprintf('>> Done!\n');
    end
    delete(gcp('nocreate'));
end

%% save the results
save(['Check_' Session],'ephys_data','ephys_data0')

%% Plot the odd and even response profiles
[~,~] = mkdir(['Check_' Session]);
cd(['Check_' Session]);

for elec=1:1024
    %fprintf(['Figure ' num2str(elec) '\n'])
    f=figure;
    hold on;
    plot(ephys_data{1}(elec,:),'r')
    plot(ephys_data{2}(elec,:),'b')
    plot(ephys_data0{1}(elec,:),'k')
    title(['Electrode ' num2str(elec)])
    l=legend({'odd','even','all'});
    set(f,'Position',[ 50 50 600 400]);
    
    saveas(f,['Elec_' num2str(elec) '.png'])
    close(f);
end
cd ..

%% Plot the odd and even response profiles === fixed ===
[~,~] = mkdir(['Check2_' Session]);
cd(['Check2_' Session]);

inv_idx = [150:-1:121 180:-1:151 210:-1:181 240:-1:211];

ephys_data2{1}=[]; ephys_data2{2}=[]; ephys_data02{1}=[];
for elec = 1:size(ephys_data{1},1)
    ephys_data2{1} = cat(1,ephys_data2{1},...
        mean([ephys_data{1}(elec,1:120);ephys_data{2}(elec,inv_idx)],1));
    ephys_data2{2} = cat(1,ephys_data2{2},...
        mean([ephys_data{2}(elec,1:120);ephys_data{1}(elec,inv_idx)],1));
    ephys_data02{1} = cat(1,ephys_data02{1},...
        mean([ephys_data0{1}(elec,1:120);ephys_data0{1}(elec,inv_idx)],1));
end  

for elec=31%:size(ephys_data2{1},1)
    %fprintf(['Figure ' num2str(elec) '\n'])
    f=figure;
    hold on;
    plot(ephys_data2{1}(elec,:),'r')
    plot(ephys_data2{2}(elec,:),'b')
    plot(ephys_data02{1}(elec,:),'k')
    title(['Electrode ' num2str(elec)])
    l=legend({'odd','even','all'});
    set(f,'Position',[ 50 50 600 400]);
    %saveas(f,['Elec_' num2str(elec) '.png'])
    %close(f);
end
cd ..

