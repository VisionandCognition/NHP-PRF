% The problem lies in the xval data where the two traces are too dissinilar
% 1) look at trace generation
%   a) Baseline/normalization
%   b) Fixation breaks

% biggest diversion even @ bar position 200-220
% also a peak in even around 170


%% Init ===================================================================
Do.Load.Any         = true;
Do.Load.ProcMUA     = true;
Do.Load.ProcLFP     = false;
Do.Load.ProcBEH     = true;
Do.SaveMUA_perArray = false;
Do.SaveLFP_perArray = false;

%% data location ==========================================================
base_path = pwd;
if ismac % laptop & portable HDD
    %data_path = '/Volumes/MRI_2/PRF_EPHYS';
    data_path = '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/Data/ephys';
    home_path = '/Users/chris';
else %
    data_path = '/media/NETDISKS/VS02/VandC/PRF_EPHYS';
    home_path = '/home/chris';
end
raw_fld = fullfile(data_path,'Data_raw');
pp_fld = fullfile(data_path,'Data_preproc');
beh_fld = fullfile(data_path,'Log_pRF');
save_fld = fullfile(data_path,'Data_check');
ch_fld = fullfile(data_path,'Channelmaps');

[~,~] = mkdir(save_fld);

%% load blackrock tools ===================================================
addpath(genpath(fullfile(home_path,'Documents','CODE','NPMK')));

%% load chronux for lfp analysis ==========================================
addpath(genpath(fullfile(home_path,'Dropbox','MATLAB_NONGIT','TOOLBOX','chronux')));

%% Load already processed files ===========================================
warning off
if Do.Load.ProcMUA
    fprintf(['Loading ' fullfile(save_fld,subj,sess,'MUA',...
        [subj '_' sess '_uncut_MUA']) '\n']);
    load(fullfile(save_fld,subj,sess,'MUA',...
        [subj '_' sess '_uncut_MUA']));
end

if Do.Load.ProcLFP
    fprintf(['Loading ' fullfile(save_fld,subj,sess,'LFP',...
        [subj '_' sess '_uncut_LFP']) '\n']);
    load(fullfile(save_fld,subj,sess,'LFP',...
        [subj '_' sess '_uncut_LFP']));
end

if Do.Load.ProcBEH
    fprintf(['Loading ' fullfile(save_fld,subj,sess,[subj '_' sess '_STIM']) '\n']);
    load(fullfile(save_fld,subj,sess,[subj '_' sess '_STIM']));
end
warning on

%% Sync timelines =========================================================
if Do.SyncTimes
    fprintf('=== Syncing time-bases across recordings and logs ===\n');
    % get log-time of all and first 'newpos'
    fprintf('Getting log-times of all and first ''newpos''\n');
    for b=1:length(B)
        log_t = [B(b).Log.Events.t];
        B(b).log_newpos = arrayfun(@(x) strcmp(x.type,'NewPosition'),B(b).Log.Events);
        B(b).t_allnp=log_t(B(b).log_newpos);
        B(b).allnp_img=[B(b).Log.Events(B(b).log_newpos).StimName];
        B(b).t_firstnp=B(b).t_allnp(1);
    end
    % get rec-times of all 'newpos'
    fprintf('Getting rec-times of all ''newpos''\n');
    for n=1:length(N)
        N(n).tsec = single(N(n).TimeStamp)./N(n).sf;
        npi = find(N(n).bits==1);
        N(n).np_samp = N(n).TimeStamp(npi);
        N(n).np_sec = N(n).tsec(npi);
        
        % split runs
        for b=1:length(B)
            ii = (b-1)*length(B(b).allnp_img)+1 : ...
                b*length(B(b).allnp_img);
            N(n).run(b).np_sec = N(n).np_sec(ii);
            N(n).run(b).barpos = B(b).allnp_img;
            N(n).run(b).stim_sec = N(n).run(b).np_sec(N(n).run(b).barpos~=0);
        end
        
    end
    fprintf('Done!\n');
end

%% Get MUA responses for each stim-position / channel =====================
if Do.SaveMUA_perArray
    win = [0 B(1).Par.TR*B(1).StimObj.Stm.RetMap.TRsPerStep]; %[start stop] in sec
    fprintf('=== Getting MUA response per trial =====\n');
    for m=1:length(M) %>> while testing, do only 1 array
        fprintf('- Array %d -', m);
        fprintf('\nChannel: 000');
        for c=1:length(M(m).chan)
            fprintf('\b\b\b\b%3d\t',c);
            M(m).ch(c).coll = [];  M(m).ch(c).coll_odd = [];  M(m).ch(c).coll_even = [];
            M(m).ch(c).BL = []; M(m).ch(c).BL_odd = []; M(m).ch(c).BL_even = [];
            
            % collect traces and average
            for r=1:length(N(m).run)
                start_samp = find(M(m).tsec >= N(m).run(r).stim_sec(1),1,'first');
                stop_samp = find(M(m).tsec >= N(m).run(r).stim_sec(end)+win(2),1,'first');
                if r==1
                    nsamp = stop_samp-start_samp;
                end
                M(m).ch(c).coll = [M(m).ch(c).coll; ...
                    M(m).chan{c}(start_samp:start_samp+nsamp)'];
                M(m).run(r).tsec = M(m).tsec(start_samp:start_samp+nsamp);
                M(m).ch(c).BL = [M(m).ch(c).BL; ...
                    M(m).chan{c}(start_samp-1000:start_samp)'];
                if mod(r,2) % odd
                    M(m).ch(c).coll_odd = [M(m).ch(c).coll_odd; ...
                        M(m).chan{c}(start_samp:start_samp+nsamp)'];
                    M(m).ch(c).BL_odd = [M(m).ch(c).BL_odd; ...
                        M(m).chan{c}(start_samp-1000:start_samp)'];
                else % even
                    M(m).ch(c).coll_even = [M(m).ch(c).coll_even; ...
                        M(m).chan{c}(start_samp:start_samp+nsamp)'];
                    M(m).ch(c).BL_even = [M(m).ch(c).BL_even; ...
                        M(m).chan{c}(start_samp-1000:start_samp)'];
                end
            end
            
            % average over runs
            mMUA(c).runs = M(m).ch(c).coll;
            mMUA(c).median = median(M(m).ch(c).coll);
            mMUA(c).BL = median(mean(M(m).ch(c).BL));
            
            mMUA_odd(c).runs = M(m).ch(c).coll_odd;
            mMUA_odd(c).median = median(M(m).ch(c).coll_odd);
            mMUA_odd(c).BL = median(mean(M(m).ch(c).BL_odd));
            
            mMUA_even(c).runs = M(m).ch(c).coll_even;
            mMUA_even(c).median = median(M(m).ch(c).coll_even);
            mMUA_even(c).BL = median(mean(M(m).ch(c).BL_even));
            
            % split by bar position
            for b=1:length(N(m).run(1).stim_sec)
                t1i= find(M(m).run(1).tsec >= ...
                    N(m).run(1).stim_sec(b)+win(1),1,'first');
                t2i= find(M(m).run(1).tsec <= ...
                    N(m).run(1).stim_sec(b)+win(2),1,'last');
                
                act_chunk = mMUA(c).median(t1i:t2i);
                mMUA(c).bar(b) = mean(act_chunk);
                clear act_chunk
                
                act_chunk_odd = mMUA_odd(c).median(t1i:t2i);
                mMUA_odd(c).bar(b) = mean(act_chunk_odd);
                clear act_chunk_odd
                
                act_chunk_even = mMUA_even(c).median(t1i:t2i);
                mMUA_even(c).bar(b) = mean(act_chunk_even);
                clear act_chunk_even
            end
        end
        fprintf('\n');
        fprintf('Saving the averaged MUA responses for array %d...\n', m);
        warning off; mkdir(fullfile(save_fld,subj,sess,'MUA')); warning on
        save(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_array_' num2str(m) '_mMUA']),'medMUA','C','-v7.3');
        clear mMUA
        save(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_array_' num2str(m) '_mMUA_odd']),'medMUA_odd','C','-v7.3');
        clear mMUA_odd
        save(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_array_' num2str(m) '_mMUA_even']),'medMUA_even','C','-v7.3');
        clear mMUA_even
    end
    clear M
    fprintf('All done!\n');
end

%% Get LFP responses for each stim-position / channel =====================
if Do.SaveLFP_perArray
    
    win = [0 B(1).Par.TR*B(1).StimObj.Stm.RetMap.TRsPerStep]; %[start stop] in sec
    
    % chronux specs
    chron.tapers = [5 9];
    chron.Fs = 500; % sampling freq
    chron.fpass = [1 150]; % [fmin fmax]
    chron.mwin = [0.500 0.050]; % moving window
    
    fprintf('=== Getting LFP response per trial =====\n');
    for m=1:length(L) %>> while testing, do only 1 array
        fprintf('- Array %d -', m);
        fprintf('\nChannel: 000');
        
        L(m).freq(1).name = 'Theta (4-8)';
        L(m).freq(1).fband = [3.9 8];
        L(m).freq(2).name = 'Alpha (8-16)';
        L(m).freq(2).fband = [8 16];
        L(m).freq(3).name = 'Beta (16-30)';
        L(m).freq(3).fband = [16 30];
        L(m).freq(4).name = 'Low gamma (30-60)';
        L(m).freq(4).fband = [30 60];
        L(m).freq(5).name = 'High gamma (60-120)';
        L(m).freq(5).fband = [60 120];
        
        for c=1:length(L(m).chan)
            fprintf('\b\b\b\b%3d\t',c);
            
            % get power traces per frequency range
            [L(m).spect(c).S,L(m).spect(c).t,L(m).spect(c).f] = ...
                mtspecgramc(L(m).chan{c}.data,chron.mwin,chron);
            L(m).spect(c).tsec = L(m).spect(c).t+L(m).t(1);
            for fb = 1:length(L(m).freq)
                fsel = L(m).spect(c).f > L(m).freq(fb).fband(1) & ...
                    L(m).spect(c).f<L(m).freq(fb).fband(2);
                L(m).ch(c).freq(fb).trace = ...
                    mean(L(m).spect(c).S(:,fsel),2);
                L(m).ch(c).freq(fb).coll = [];
                L(m).ch(c).freq(fb).BL = [];
                L(m).ch(c).freq(fb).coll_odd = [];
                L(m).ch(c).freq(fb).BL_odd = [];
                L(m).ch(c).freq(fb).coll_even = [];
                L(m).ch(c).freq(fb).BL_even = [];
            end
            
            % collect traces and average
            for r=1:length(N(m).run)
                start_samp = find(L(m).spect(c).tsec >= ...
                    N(m).run(r).stim_sec(1),1,'first');
                stop_samp = find(L(m).spect(c).tsec >= ...
                    N(m).run(r).stim_sec(end)+win(2),1,'first');
                
                if r==1
                    nsamp = stop_samp-start_samp;
                end
                
                for fb = 1:length(L(m).freq)
                    L(m).ch(c).freq(fb).coll = [L(m).ch(c).freq(fb).coll; ...
                        L(m).ch(c).freq(fb).trace(start_samp:start_samp+nsamp)'];
                    L(m).ch(c).freq(fb).BL = [L(m).ch(c).freq(fb).BL; ...
                        L(m).ch(c).freq(fb).trace(start_samp-500:start_samp)'];
                    
                    if mod(r,2) % odd
                        L(m).ch(c).freq(fb).coll_odd = [L(m).ch(c).freq(fb).coll_odd; ...
                            L(m).ch(c).freq(fb).trace(start_samp:start_samp+nsamp)'];
                        L(m).ch(c).freq(fb).BL_odd = [L(m).ch(c).freq(fb).BL_odd; ...
                            L(m).ch(c).freq(fb).trace(start_samp-500:start_samp)'];
                    else % even
                        L(m).ch(c).freq(fb).coll_even = [L(m).ch(c).freq(fb).coll_even; ...
                            L(m).ch(c).freq(fb).trace(start_samp:start_samp+nsamp)'];
                        L(m).ch(c).freq(fb).BL_even = [L(m).ch(c).freq(fb).BL_even; ...
                            L(m).ch(c).freq(fb).trace(start_samp-500:start_samp)'];
                        
                    end
                    
                    
                end
                L(m).run(r).tsec = L(m).spect(c).tsec(start_samp:start_samp+nsamp);
            end
            
            % average over runs
            for fb = 1:length(L(m).freq)
                mLFP(c).freq(fb).runs = L(m).ch(c).freq(fb).coll;
                mLFP(c).freq(fb).mean = mean(L(m).ch(c).freq(fb).coll);
                mLFP(c).freq(fb).std = std(L(m).ch(c).freq(fb).coll);
                mLFP(c).freq(fb).BL = mean(mean(L(m).ch(c).freq(fb).BL));
                
                mLFP_odd(c).freq(fb).runs = L(m).ch(c).freq(fb).coll_odd;
                mLFP_odd(c).freq(fb).mean = mean(L(m).ch(c).freq(fb).coll_odd);
                mLFP_odd(c).freq(fb).std = std(L(m).ch(c).freq(fb).coll_odd);
                mLFP_odd(c).freq(fb).BL = mean(mean(L(m).ch(c).freq(fb).BL_odd));
                
                mLFP_even(c).freq(fb).runs = L(m).ch(c).freq(fb).coll_even;
                mLFP_even(c).freq(fb).mean = mean(L(m).ch(c).freq(fb).coll_even);
                mLFP_even(c).freq(fb).std = std(L(m).ch(c).freq(fb).coll_even);
                mLFP_even(c).freq(fb).BL = mean(mean(L(m).ch(c).freq(fb).BL_even));
                
                % split by bar position
                for b=1:length(N(m).run(1).stim_sec)
                    t1i= find(L(m).run(1).tsec >= ...
                        N(m).run(1).stim_sec(b)+win(1),1,'first');
                    t2i= find(L(m).run(1).tsec <= ...
                        N(m).run(1).stim_sec(b)+win(2),1,'last');
                    
                    act_chunk = mLFP(c).freq(fb).mean(t1i:t2i);
                    mLFP(c).freq(fb).bar(b) = mean(act_chunk);
                    clear act_chunk
                    
                    act_chunk_odd = mLFP_odd(c).freq(fb).mean(t1i:t2i);
                    mLFP_odd(c).freq(fb).bar(b) = mean(act_chunk_odd);
                    clear act_chunk_odd
                    
                    act_chunk_even = mLFP_even(c).freq(fb).mean(t1i:t2i);
                    mLFP_even(c).freq(fb).bar(b) = mean(act_chunk_even);
                    clear act_chunk_even
                end
                
            end
            
        end
        fprintf('\n');
        fprintf('Saving the averaged LFP responses for array %d...\n', m);
        warning off; mkdir(fullfile(save_fld,subj,sess,'LFP')); warning on
        save(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_array_' num2str(m) '_mLFP']),'mLFP','C','-v7.3');
        clear mLFP
        warning off; mkdir(fullfile(save_fld,subj,sess,'LFP')); warning on
        save(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_array_' num2str(m) '_mLFP_odd']),'mLFP_odd','C','-v7.3');
        clear mLFP_odd
        warning off; mkdir(fullfile(save_fld,subj,sess,'LFP')); warning on
        save(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_array_' num2str(m) '_mLFP_even']),'mLFP_even','C','-v7.3');
        clear mLFP_even
    end
    clear L
    fprintf('All done!\n');
end