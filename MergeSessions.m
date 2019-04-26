%% Merge median BOLD signals
clear all; clc;
MONKEY={'danny','eddy'};
for m=1:length(MONKEY)
    fprintf(['==== Processing monkey: ' MONKEY{m} ' ====\n']);
    cd(['pRF_sub-' MONKEY{m} '_us-padded']);
    fls = dir('medianBOLD*');
    MB.medianBOLD=[]; MB.medianBOLD_inv=[];
    MB.stim={}; MB.stim_inv={};
    for f=1:length(fls)
        fprintf(['Adding ' fls(f).name '\n']);
        M = load(fls(f).name);
        MB(f).medianBOLD = M.medianBOLD;
        MB(f).stim = M.stim;
        if isfield(M,'medianBOLD_inv')
            MB(f).medianBOLD_inv = M.medianBOLD_inv;
            MB(f).stim_inv = M.stim_inv;
        end
        clear M
    end
    save('AllSessions','MB');
    cd ..
    clear MB
    fprintf('\n\n');
end


%% Average (or median)
