function ck_GetEphys_pRF2(monkeys,models,output,dataset)

if nargin < 2
    fprintf('ERROR: Not enough arguments specified\n');
    return
end
clc;

%% data location ==========================================================
fitres_path = ...
    ['/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/' ...
    'pRF/FitResults/ephys/Results_' dataset];
save_path = ...
    ['/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/' ...
    'pRF/FitResults/ephys/Combined'];
chanmap_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/Data/ephys';
start_fld=pwd;

%% Collect EPHYS ==========================================================
fprintf('Collecting ephys retinotopic maps...\n');

for m = 1:length(monkeys)
    fprintf(['Processing monkey: ' monkeys{m} '\n']);
    cMonkey = monkeys{m}; cMonkey(1)=upper(cMonkey(1));
    
    R(m).monkey = monkeys{m};
    R(m).mode = 'ephys';
    
    fprintf(['Getting channel maps \n']);
    CM = load(fullfile(chanmap_path,monkeys{m},...
        ['channel_area_mapping_' monkeys{m}]));
    R(m).ChanMap = CM; clear CM;
    
    for mm = 1:length(models)
        fprintf(['Model: ' models{mm} '\n']);
        R(m).model(mm).prfmodel = models{mm};
        
        if ~strcmp(models{mm},'classicRF')
            for i = 1:8
                fprintf(['Instance: ' num2str(i) '\n']);
                %% MUA --
                fprintf('MUA\n');
                RES = load(fullfile(fitres_path,monkeys{m},models{mm},...
                    ['Instance_' num2str(i)],...
                    ['pRF_Sess-mua_Inst_' num2str(i)]),'result');
                R(m).model(mm).MUA(i) = RES.result;
                clear RES
                
                %% LFP --
                for l=1:5
                    fprintf(['LFP ' num2str(l) '\n']);
                    RES = load(fullfile(fitres_path,monkeys{m},models{mm},...
                        ['Instance_' num2str(i)],...
                        ['pRF_Sess-lfp_fb' num2str(l) '_Inst_' num2str(i)]),'result');
                    R(m).model(mm).LFP(i,l) = RES.result;
                    clear RES
                end
            end
        else
            for i = 1:8
                fprintf(['Instance: ' num2str(i) '\n']);
                RES = load(fullfile(fitres_path,monkeys{m},models{mm},...
                    ['RFs_instance' num2str(i)]));
                R(m).model(mm).MUA(i).RF = RES.RFs;
                R(m).model(mm).MUA(i).chanRF = RES.channelRFs;
                R(m).model(mm).MUA(i).SNR = RES.meanChannelSNR;
                clear RES
            end
        end
    end  
end

%% Save the combined results ==============================================
fprintf('Saving the combined ephys result-file\n');
[~,~,~] = mkdir(fullfile(save_path,dataset));
save(fullfile(save_path,dataset,output),'R')
