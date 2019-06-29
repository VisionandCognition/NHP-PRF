%% Merge median BOLD signals
clear all; clc;
MONKEY={'danny','eddy'};
startfld=pwd;
for m=1:length(MONKEY)
    fprintf(['==== Processing monkey: ' MONKEY{m} ' ====\n']);
    cd ..; cd Data;
    cd(['pRF_sub-' MONKEY{m} '_us-padded']);
    fls = dir('medianBOLD*');
    
    MB.medianBOLD=[]; MB.medianBOLD_inv=[];
    MB.nRuns=[]; MB.nRuns_inv=[];
    MB.stim={}; MB.stim_inv={};
    
    for f=1:length(fls)
        fprintf(['Adding ' fls(f).name '\n']);
        M = load(fls(f).name);
        MB.medianBOLD = cat(5,MB.medianBOLD,M.medianBOLD);
        MB.nRuns = cat(5,MB.nRuns,M.nRuns);
        MB.stim = M.stim;
        if isfield(M,'medianBOLD_inv')
            MB.medianBOLD_inv = cat(5,MB.medianBOLD_inv,M.medianBOLD_inv);
            MB.nRuns_inv = cat(5,MB.nRuns_inv,M.nRuns_inv);
            MB.stim_inv = M.stim_inv;
        end
        clear M
    end
    stim.norm=MB.stim;
    % normal mean
    sess_meanBOLD = nanmean(MB.medianBOLD,5);
    % weighted mean by number of runs
    sess_wmeanBOLD = nansum(MB.medianBOLD.*MB.nRuns,5)./...
        nansum(MB.nRuns,5);
    % std
    sess_sdBOLD = nanstd(MB.medianBOLD,1,5);
    %median
    sess_medianBOLD = nanmedian(MB.medianBOLD,5);
    if isfield(MB,'medianBOLD_inv')
        stim.inv=MB.stim_inv;
        % normal mean
        sess_meanBOLD_inv = nanmean(MB.medianBOLD_inv,5);
        % weighted mean by number of runs
        sess_wmeanBOLD_inv = nansum(MB.medianBOLD_inv.*MB.nRuns_inv,5)./...
            nansum(MB.nRuns_inv,5);
        % std
        sess_sdBOLD_inv = nanstd(MB.medianBOLD_inv,1,5);
        %median
        sess_medianBOLD_inv = nanmedian(MB.medianBOLD_inv,5);
    end
    save('AllSessions-avg','stim','sess_meanBOLD','sess_meanBOLD_inv',...
        'sess_wmeanBOLD','sess_wmeanBOLD_inv','sess_medianBOLD','sess_medianBOLD_inv',...
        'sess_sdBOLD','sess_sdBOLD_inv','MB');
    save('AllSessions-only_avg','stim','sess_meanBOLD','sess_meanBOLD_inv',...
        'sess_wmeanBOLD','sess_wmeanBOLD_inv','sess_medianBOLD','sess_medianBOLD_inv',...
        'sess_sdBOLD','sess_sdBOLD_inv');
    cd(startfld);
    clear MB
    fprintf('\n\n');
end