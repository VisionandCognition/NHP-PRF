function pRF_avg_BOLD_tseries(monkey,sess)

%% load the equal length files
% manual for now
cd(['pRF_sub-' monkey '_us-padded']);
load(['ses-' sess '-230vols']); %#ok<*LOAD>
fprintf(['Processing ses-' sess '-230vols\n']); %#ok<*LOAD>

%% average
% stim normal of inverted?
r_fw=[]; r_inv=[];
for r=1:length(p_run)
    if p_run(r).stim{11}(10,80)
        r_fw=[r_fw r];
    else
        r_inv=[r_inv r];
    end
    
%     if ~isempty(p_run(r).stimori)
%         so=s_run(r).stimori(6,2);
%     else
%         so=s_run(1).stimori(6,2);
%     end
%     if  so == 90
%         r_fw=[r_fw r]; %#ok<*AGROW>
%     elseif so == 270
%         r_inv=[r_inv r];
%     else
%         fprintf('ERROR: unknown stimulus configuration');
%     end
end

%% FW stim
fprintf('Adressing regular stimuli \n')
% pre-allocate collection arrays
run_fw=nan([size(p_run(1).vol{1}) length(p_run(r).vol) length(r_fw)]);
% collect
for r=r_fw
    fprintf(['Processing r = ' num2str(r) '\n']);
    nanvol = nan(size(p_run(r).vol{1}));
    vol = p_run(r).vol;
    %vol(p_run(r).excvol) = {nanvol};
    for t=1:length(vol)
        run_fw(:,:,:,t,r) = vol{t};
    end
    % convert to percentage BOLD change
    mSig=nanmean(run_fw(:,:,:,:,r),4);
    for d4=1:size(run_fw,4)
        NormVol=100.*((run_fw(:,:,:,d4,r)-mSig)./mSig);
        NormVol(isnan(NormVol))=0;
        run_fw(:,:,:,d4,r)=NormVol;
    end
    % select timepoints to include
    for v=find(p_run(r).inc==0)
        run_fw(:,:,:,v,r)=nanvol;
    end
end
fprintf('Getting the median BOLD signal for all voxels\n');
medianBOLD = nanmedian(run_fw,5); %#ok<*NASGU>
stim = p_run(r).stim;
clear run_fw

%% INV stim
fprintf('Adressing inverse stimuli\n')
if ~isempty(r_inv)
    % pre-allocate collection arrays
    run_inv=nan([size(p_run(1).vol{1}) length(p_run(r).vol) length(r_inv)]);
    % collect
    for r=r_inv
        fprintf(['Processing r = ' num2str(r) '\n']);
        nanvol = nan(size(p_run(r).vol{1}));
        vol = p_run(r).vol;
        %vol(p_run(r).excvol) = {nanvol};
        for t=1:length(vol)
            run_inv(:,:,:,t,r) = vol{t};
        end
        % convert to percentage BOLD change
        mSig=nanmean(run_inv(:,:,:,:,r),4);
        for d4=1:size(run_inv,4)
            NormVol=100.*((run_inv(:,:,:,d4,r)-mSig)./mSig);
            NormVol(isnan(NormVol))=0;
            run_inv(:,:,:,d4,r)=NormVol;
        end
        % select timepoints to include
        for v=find(p_run(r).inc==0)
            run_inv(:,:,:,v,r)=nanvol;
        end
    end
    fprintf('Getting the median BOLD signal for all voxels\n');
    medianBOLD_inv = nanmedian(run_inv,5);
    stim_inv = p_run(r).stim;
    clear run_inv
else
    fprintf('Does not exist for this session\n')
end

%% save
fprintf('Saving the result\n');
if ~isempty(r_inv)
    save(['medianBOLD_sess-' sess],'medianBOLD','medianBOLD_inv','stim','stim_inv','-v7.3');
else
    save(['medianBOLD_sess-' sess],'medianBOLD','stim','-v7.3');
end
cd ..