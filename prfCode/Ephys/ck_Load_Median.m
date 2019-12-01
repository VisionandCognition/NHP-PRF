function ck_Load_Median(subj, sess, Do)
% loads the (Blackroack) ephys PRF data and organizes it in trials
% c.klink@nin.knaw.nl

%% init ===================================================================
if nargin<3
    fprintf('Not specified what to do. Taking defaults\n')
    Do.Load.Any         = true;
    Do.Load.DigChan     = false; % new files
    Do.Load.MUA         = false; % new files
    Do.Load.LFP         = false; % new files
    Do.Load.Behavior    = false; % new files
    
    Do.Load.ProcMUA     = true;
    Do.Load.ProcLFP     = true;
    Do.Load.ProcBEH     = true;
    
    Do.SyncTimes        = true;
    Do.SaveUncut        = false;
    
    Do.SaveMUA_perArray = true;
    Do.SaveLFP_perArray = true;
    
    Do.CreatePrediction = false;
    
    Do.FitPRF           = false;
    Do.FitMUA           = false; % MUA data
    Do.FitLFP           = false; % LFP data (split by freq band)
    
    Do.PlotPRF_MUA      = false;
    Do.PlotPRF_LFP      = false;
end

if nargin<2
    subj = 'Lick';
    sess = '20180807_B2';
    fprintf(['Using defaults, SUBJ: ' subj ', SESS: ' sess '\n\n']);
end
setup = 'NIN';

fprintf('=============================================================\n');
fprintf(' Loading and processing data for %s, %s\n', subj, sess);
fprintf('=============================================================\n');

%% data location ==========================================================
base_path = pwd;
if ismac % laptop & portable HDD
    data_path = '/Volumes/MRI_2/PRF_EPHYS';
    home_path = '/Users/chris';
else %
    data_path = '/media/NETDISKS/VS02/VandC/PRF_EPHYS';
    home_path = '/home/chris';
end
raw_fld = fullfile(data_path,'Data_raw');
pp_fld = fullfile(data_path,'Data_preproc');
beh_fld = fullfile(data_path,'Log_pRF');
save_fld = fullfile(data_path,'Data_proc');
ch_fld = fullfile(data_path,'Channelmaps');

%% load blackrock tools ===================================================
addpath(genpath(fullfile(home_path,'Documents','CODE','NPMK')));

%% load chronux for lfp analysis ==========================================
addpath(genpath(fullfile(home_path,'Dropbox','MATLAB_NONGIT','TOOLBOX','chronux')));

%% Get digital channel contents ===========================================
if Do.Load.DigChan
    cd(fullfile(raw_fld,subj,sess));
    nev_files = dir('*.nev');
    fprintf('=== Loading digital channel info ===\n');
    for i=1:length(nev_files)
        NEV=openNEV(fullfile(nev_files(i).folder,nev_files(i).name));
        N(i).TimeStamp = NEV.Data.SerialDigitalIO.TimeStamp;
        N(i).t0 = N(i).TimeStamp(1);
        
        N(i).DigData = NEV.Data.SerialDigitalIO.UnparsedData;
        N(i).bits = log2(single(N(i).DigData));
        
        N(i).sf = 30e3;
        N(i).fp_i = find(N(i).bits==1,1,'first');
        N(i).t_fp = N(i).TimeStamp(N(i).fp_i);
        clear NEV
    end
    fprintf('Done!\n');
    cd(base_path);
end

% get the channelmapping
if strcmp(subj,'Lick')
    C=load(fullfile(ch_fld,'channel_area_mapping_lick'));
elseif strcmp(subj,'Aston')
    C=load(fullfile(ch_fld,'channel_area_mapping_aston'));
end

%% Get MUA ================================================================
if Do.Load.MUA
    cd(fullfile(pp_fld,subj,sess));
    mua_files = dir('MUA*.mat');
    fprintf('=== Loading preprocessed MUA ===\n');
    for i=1:length(mua_files)
        fprintf([num2str(i) '/' num2str(length(mua_files)) '\n']);
        load(fullfile(mua_files(i).folder,mua_files(i).name));
        M(i).chan=channelDataMUA;% (1:128);
        M(i).sf = 1e3;
        % create timeline based on digital channel
        M(i).tsec = 0:1/M(i).sf:(1/M(i).sf)*length(M(i).chan{1});
        M(i).tsec(end)=[];
        clear channelDataMUA
    end
    fprintf('Done!\n');
    cd(base_path);
end

%% Get LFP ================================================================
if Do.Load.LFP
    cd(fullfile(pp_fld,subj,sess));
    lfp_files = dir('LFP*.mat');
    fprintf('=== Loading preprocessed LFP ===\n');
    for i=1:length(lfp_files)
        fprintf([num2str(i) '/' num2str(length(lfp_files)) '\n']);
        load(fullfile(lfp_files(i).folder,lfp_files(i).name));
        L(i).chan=channelDataLFP(1:128);
        L(i).t = M(i).tsec(1):1/L(i).chan{1}.LFPsamplerate:M(i).tsec(end);
        clear channelDataLFP
    end
    fprintf('Done!\n');
    cd(base_path);
end

%% Get behavioral log =====================================================
if Do.Load.Behavior
    cd(fullfile(beh_fld,setup,subj))
    temp_fld=dir(['*' sess(1:8) '*']);
    cd(fullfile(temp_fld(1).folder,temp_fld(1).name));
    temp_fld=dir(['*' sess(1:8) '*']);
    fprintf('=== Loading Behavioral logs ===\n');
    for recn = 1:length(temp_fld)
        cd(fullfile(temp_fld(recn).folder,temp_fld(recn).name));
        mat_fn = dir('Log*Run*.mat');
        for m=1:length(mat_fn)
            fprintf([mat_fn(m).name '\n']);
            B(recn,m)=load(fullfile(mat_fn(m).folder,mat_fn(m).name));
            MTi = find(arrayfun(@(x) strcmp(x.type,'MRITrigger'),...
                B(recn,m).Log.Events)==1,1,'first');
            if m==1
                triggertime=B(recn,m).Log.Events(MTi).t;
            end
            for tt=1:length(B(recn,m).Log.Events)
                B(recn,m).Log.Events(tt).t_mri = B(recn,m).Log.Events(tt).t - ...
                    triggertime;
            end
        end
        % also get the stimulus
        S = load('RetMap_Stimulus.mat');
        S = S.ret_vid;
        
        for j=1:size(S,2)
            ori = S(j).orient(1);
            imnr = mod(j,B(recn).StimObj.Stm.RetMap.nSteps/8);
            if imnr==0; imnr=B(recn).StimObj.Stm.RetMap.nSteps/8;end
            img = S(imnr).img{1};
            img(img==255 | img==0)=1;
            img(img~=1)=0;
            STIM.img{j}=imrotate(img,-ori,'nearest','crop');
        end
        STIM.blank = uint8(zeros(size(STIM.img{j})));
        fprintf(['Saving ' fullfile(save_fld,[subj '_' sess '_STIM\n'])]);
        warning off; mkdir(fullfile(save_fld,subj,sess)); warning on
        save(fullfile(save_fld,subj,sess,[subj '_' sess '_STIM']),'STIM','B','-v7.3');
    end
    fprintf('Done!\n');
    cd(base_path);
end

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

%% Saving the (rawish) MUA & LFP data =====================================
if Do.SaveUncut
    fprintf('=== Saving the un-cut data =====\n');
    fprintf(fullfile(save_fld,subj,sess,'MUA',...
        [subj '_' sess '_uncut_MUA\n']));
    warning off; mkdir(fullfile(save_fld,subj,sess,'MUA')); warning on
    save(fullfile(save_fld,subj,sess,'MUA',...
        [subj '_' sess '_uncut_MUA']),'M','B','N','C','-v7.3');
    fprintf(fullfile(save_fld,subj,sess,'LFP',...
        [subj '_' sess '_uncut_LFP\n']));
    warning off; mkdir(fullfile(save_fld,subj,sess,'LFP')); warning on
    save(fullfile(save_fld,subj,sess,'LFP',...
        [subj '_' sess '_uncut_LFP']),'L','B','N','C','-v7.3');
    fprintf('Done!\n');
end

%% Get MUA responses for each stim-position / channel =====================
if Do.SaveMUA_perArray
    win = [0.05 B(1).Par.TR*B(1).StimObj.Stm.RetMap.TRsPerStep]; %[start stop] in sec
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
            
%             % ----- OLD VERSION -------------------------------------------
%             % average over runs
%             mMUA(c).runs = M(m).ch(c).coll;
%             mMUA(c).mean = mean(M(m).ch(c).coll);
%             mMUA(c).std = std(M(m).ch(c).coll);
%             mMUA(c).BL = mean(mean(M(m).ch(c).BL));
%             
%             mMUA_odd(c).runs = M(m).ch(c).coll_odd;
%             mMUA_odd(c).mean = mean(M(m).ch(c).coll_odd);
%             mMUA_odd(c).std = std(M(m).ch(c).coll_odd);
%             mMUA_odd(c).BL = mean(mean(M(m).ch(c).BL_odd));
%             
%             mMUA_even(c).runs = M(m).ch(c).coll_even;
%             mMUA_even(c).mean = mean(M(m).ch(c).coll_even);
%             mMUA_even(c).std = std(M(m).ch(c).coll_even);
%             mMUA_even(c).BL = mean(mean(M(m).ch(c).BL_even));
%             
%             % split by bar position
%             for b=1:length(N(m).run(1).stim_sec)
%                 t1i= find(M(m).run(1).tsec >= ...
%                     N(m).run(1).stim_sec(b)+win(1),1,'first');
%                 t2i= find(M(m).run(1).tsec <= ...
%                     N(m).run(1).stim_sec(b)+win(2),1,'last');
%                 
%                 act_chunk = mMUA(c).mean(t1i:t2i);
%                 mMUA(c).bar(b) = mean(act_chunk);
%                 clear act_chunk
%                 
%                 act_chunk_odd = mMUA_odd(c).mean(t1i:t2i);
%                 mMUA_odd(c).bar(b) = mean(act_chunk_odd);
%                 clear act_chunk_odd
%                 
%                 act_chunk_even = mMUA_even(c).mean(t1i:t2i);
%                 mMUA_even(c).bar(b) = mean(act_chunk_even);
%                 clear act_chunk_even
%             end
%             % -------------------------------------------------------------
            
            % median over runs
            mMUA(c).runs = M(m).ch(c).coll;
            SIG = [];
            for r=1:size(M(m).ch(c).BL,1)
                 SIG = [SIG; M(m).ch(c).coll(r,:)-...
                     mean(M(m).ch(c).BL(r,:),2)];
            end
            mMUA(c).median = median(SIG);
            mMUA(c).BL = median(M(m).ch(c).BL);
            
            mMUA_odd(c).runs = M(m).ch(c).coll_odd;
            SIG_odd = [];
            for r=1:size(M(m).ch(c).BL_odd,1)
                 SIG_odd = [SIG_odd; M(m).ch(c).coll_odd(r,:)-...
                     mean(M(m).ch(c).BL_odd(r,:),2)];
            end
            mMUA_odd(c).median = median(SIG_odd);
            mMUA_odd(c).BL = median(M(m).ch(c).BL_odd);
            
            mMUA_even(c).runs = M(m).ch(c).coll_even;
            SIG_even = [];
            for r=1:size(M(m).ch(c).BL_even,1)
                 SIG_even = [SIG_even; M(m).ch(c).coll_even(r,:)-...
                     mean(M(m).ch(c).BL_even(r,:),2)];
            end
            mMUA_even(c).median = median(SIG_even);
            mMUA_even(c).BL = median(M(m).ch(c).BL_even);
            
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
        fprintf('Saving the median MUA responses for array %d...\n', m);
        warning off; [~,~]=mkdir(fullfile(save_fld,subj,sess,'MUA')); warning on
        save(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_array_' num2str(m) '_medMUA']),'mMUA','C','-v7.3');
        clear mMUA
        save(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_array_' num2str(m) '_medMUA_odd']),'mMUA_odd','C','-v7.3');
        clear mMUA_odd
        save(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_array_' num2str(m) '_medMUA_even']),'mMUA_even','C','-v7.3');
        clear mMUA_even
    end
    clear M
    fprintf('All done!\n');
end

%% Get LFP responses for each stim-position / channel =====================
if Do.SaveLFP_perArray
    
    win = [0.05 B(1).Par.TR*B(1).StimObj.Stm.RetMap.TRsPerStep]; %[start stop] in sec
    
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
                mLFP(c).freq(fb).median = median(L(m).ch(c).freq(fb).coll);
                mLFP(c).freq(fb).BL = median(mean(L(m).ch(c).freq(fb).BL));
                
                mLFP_odd(c).freq(fb).runs = L(m).ch(c).freq(fb).coll_odd;
                mLFP_odd(c).freq(fb).mean = median(L(m).ch(c).freq(fb).coll_odd);
                mLFP_odd(c).freq(fb).BL = median(mean(L(m).ch(c).freq(fb).BL_odd));
                
                mLFP_even(c).freq(fb).runs = L(m).ch(c).freq(fb).coll_even;
                mLFP_even(c).freq(fb).mean = median(L(m).ch(c).freq(fb).coll_even);
                mLFP_even(c).freq(fb).BL = median(mean(L(m).ch(c).freq(fb).BL_even));
                
                % split by bar position
                for b=1:length(N(m).run(1).stim_sec)
                    t1i= find(L(m).run(1).tsec >= ...
                        N(m).run(1).stim_sec(b)+win(1),1,'first');
                    t2i= find(L(m).run(1).tsec <= ...
                        N(m).run(1).stim_sec(b)+win(2),1,'last');
                    
                    act_chunk = mLFP(c).freq(fb).median(t1i:t2i);
                    mLFP(c).freq(fb).bar(b) = mean(act_chunk);
                    clear act_chunk
                    
                    act_chunk_odd = mLFP_odd(c).freq(fb).median(t1i:t2i);
                    mLFP_odd(c).freq(fb).bar(b) = mean(act_chunk_odd);
                    clear act_chunk_odd
                    
                    act_chunk_even = mLFP_even(c).freq(fb).median(t1i:t2i);
                    mLFP_even(c).freq(fb).bar(b) = mean(act_chunk_even);
                    clear act_chunk_even
                end
                
            end
            
        end
        fprintf('\n');
        fprintf('Saving the median LFP responses for array %d...\n', m);
        warning off; mkdir(fullfile(save_fld,subj,sess,'LFP')); warning on
        save(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_array_' num2str(m) '_medLFP']),'mLFP','C','-v7.3');
        clear mLFP
        warning off; mkdir(fullfile(save_fld,subj,sess,'LFP')); warning on
        save(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_array_' num2str(m) '_medLFP_odd']),'mLFP_odd','C','-v7.3');
        clear mLFP_odd
        warning off; mkdir(fullfile(save_fld,subj,sess,'LFP')); warning on
        save(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_array_' num2str(m) '_medLFP_even']),'mLFP_even','C','-v7.3');
        clear mLFP_even
    end
    clear L
    fprintf('All done!\n');
end