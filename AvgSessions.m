%% Merge median BOLD signals
clear all; clc;
MONKEY={'danny','eddy'};
for m=1%:length(MONKEY)
    fprintf(['==== Processing monkey: ' MONKEY{m} ' ====\n']);
    cd(['pRF_sub-' MONKEY{m} '_us-padded']);
    fls = dir('medianBOLD*');
    
    MB.medianBOLD=[]; MB.medianBOLD_inv=[];
    MB.stim={}; MB.stim_inv={};
    
    for f=1:length(fls)
        fprintf(['Adding ' fls(f).name '\n']);
        M = load(fls(f).name);
        MB.medianBOLD = cat(5,MB.medianBOLD,M.medianBOLD);
        MB.stim = M.stim;
        if isfield(M,'medianBOLD_inv')
            MB.medianBOLD_inv = cat(5,MB.medianBOLD_inv,M.medianBOLD_inv);
            MB.stim_inv = M.stim_inv;
        end
        clear M
    end
    stim.norm=MB.stim;
    sess_meanBOLD = nanmean(MB.medianBOLD,5);
    sess_sdBOLD = nanstd(MB.medianBOLD,1,5);
    sess_medianBOLD = nanmedian(MB.medianBOLD,5);
    if isfield(MB,'medianBOLD_inv')
        stim.inv=MB.stim_inv;
        sess_meanBOLD_inv = nanmean(MB.medianBOLD_inv,5);
        sess_sdBOLD_inv = nanstd(MB.medianBOLD_inv,1,5);
        sess_medianBOLD_inv = nanmedian(MB.medianBOLD_inv,5);
    end
    save('AllSessions-avg','stim','sess_meanBOLD','sess_meanBOLD_inv',...
        'sess_medianBOLD','sess_medianBOLD_inv',...
        'sess_sdBOLD','sess_sdBOLD_inv','MB');
    cd ..
    clear MB
    fprintf('\n\n');
end

%%
load('AllSessions-avg','stim','sess_meanBOLD','sess_meanBOLD_inv',...
        'sess_medianBOLD','sess_medianBOLD_inv',...
        'sess_sdBOLD','sess_sdBOLD_inv')
save('AllSessions-only_avg','stim','sess_meanBOLD','sess_meanBOLD_inv',...
        'sess_medianBOLD','sess_medianBOLD_inv',...
        'sess_sdBOLD','sess_sdBOLD_inv');