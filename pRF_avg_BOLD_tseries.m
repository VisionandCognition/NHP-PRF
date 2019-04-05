function pRF_avg_BOLD_tseries(monkey,sess)

%% load the equal length files
% manual for now
cd(['pRF_sub-' monkey '_us-avg']);
load(['ses-' sess '-equal_length']); %#ok<*LOAD>

%% average
% stim normal of inverted?
r_fw=[]; r_inv=[];
for r=1:length(s_run)
    if ~isempty(s_run(r).stimori)
        so=s_run(r).stimori(6,2);
    else
        so=s_run(1).stimori(6,2);
    end
    if  so == 90
        r_fw=[r_fw r]; %#ok<*AGROW>
    elseif so == 270
        r_inv=[r_inv r];
    else
        fprintf('ERROR: unknown stimulus configuration');
    end
end

%% FW stim
fprintf('Adressing regular stimuli \n')
% pre-allocate collection arrays
run_fw=nan([size(s_run(1).vol{1}) length(s_run(r).vol) length(r_fw)]);
% collect
for r=r_fw
    fprintf(['Processing r = ' num2str(r) '\n']);
    nanvol = nan(size(s_run(r).vol{1}));
    vol = s_run(r).vol;
    vol(s_run(r).excvol) = {nanvol};
    for t=1:length(vol)
        run_fw(:,:,:,t,r) = vol{t};
    end
end
fprintf('Getting the median BOLD signal for all voxels\n');
medianBOLD = nanmedian(run_fw,5); %#ok<*NASGU>
stim = s_run(r).stim;
clear run_fw

%% INV stim
fprintf('Adressing inverse stimuli\n')
if ~isempty(r_inv)
    % pre-allocate collection arrays
    run_inv=nan([size(s_run(1).vol{1}) length(s_run(r).vol) length(r_inv)]);
    % collect
    for r=r_inv
        fprintf(['Processing r = ' num2str(r) '\n']);
        nanvol = nan(size(s_run(r).vol{1}));
        vol = s_run(r).vol;
        vol(s_run(r).excvol) = {nanvol};
        for t=1:length(vol)
            run_inv(:,:,:,t,r) = vol{t};
        end
    end
    fprintf('Getting the median BOLD signal for all voxels\n');
    medianBOLD_inv = nanmedian(run_inv,5);
    stim_inv = s_run(r).stim;
    clear run_inv
else
    fprintf('Does not exist for this session\n')
end

%% save
if ~isempty(r_inv)
    save(['medianBOLD_sess-' sess],'medianBOLD','medianBOLD_inv','stim','stim_inv','-v7.3');
else
    save(['medianBOLD_sess-' sess],'medianBOLD','stim','-v7.3');
end
cd ..