function ck_FitPRF(subj, sess, Do)
% Fit the PRF prediction to the ephys data
% c.klink@nin.knaw.nl

%% init ===================================================================
if nargin<3
    Do.FitMUA = true; % MUA data
    Do.FitLFP = true; % LFP data (split by freq band)
    Do.AllowNegBetas = true; % If we allow negative beta weights we can also fit inhibition
end

if nargin<2
    subj = 'Lick';
    sess = '20180807_B2';
    fprintf(['Using defaults, SUBJ: ' subj ', SESS: ' sess '\n\n']);
end

fprintf('=============================================================\n');
fprintf(' Loading and processing data for %s, %s\n', subj, sess);
fprintf('=============================================================\n');

%% data location ==========================================================
base_path = pwd;
if ismac % laptop & portable HDD
    data_path = '/Volumes/MRI_2/PRF_EPHYS';
    home_path = '/Users/chris';
    ncores = 4; % parallel on how many cores?
elseif isunix %
    data_path = '/media/NETDISKS/VCNIN/PRF_EPHYS';
    home_path = '/home/chris';
    ncores = 6; % do this thing in parallel on how many cores?
else % windows >> running on server
    data_path = '\\vcnin\PRF_EPHYS';
    home_path = 'D:\CK\code\MATLAB\CURRENT\NIN\PHYSIOLOGY\pRF';
    ncores = 12; % do this thing in parallel on how many cores?
end
raw_fld = fullfile(data_path,'Data_raw');
pp_fld = fullfile(data_path,'Data_preproc');
beh_fld = fullfile(data_path,'Log_pRF');
save_fld = fullfile(data_path,'Data_proc');

%% Load predictions =======================================================
fprintf('Loading pRF prediction...\n');
load(fullfile(data_path,'PRF_Prediction.mat'));

%% Load MUA responses & find best fitting pRFs ============================
if Do.FitMUA
    fprintf('Loading and processing MUA pRF results...\n');
    %% cylce over arrays/instances
    for a=1:8
        %% start
        fprintf('Array %d\n', a);
        load(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_array_' num2str(a) '_mMUA']));
        
        MUAm=[];
        for c=1:length(mMUA)
            MUAm=[MUAm mMUA(c).bar'-mMUA(c).BL]; % subtract baseline
            %MUAm=[MUAm mMUA(a,c).mBarResp(5:end)]; 
        end
        
        nchans = size(mMUA,2);
        for s = 1:12
            slice(s).ind = ones(nchans,1);
            slice(s).bestbeta = NaN(2,nchans);
            slice(s).currentbest = inf(1,nchans);
        end
        
        %% Divide the timecourse up into slices
        ngauss = size(PRF_pred,2);
        nslices = 12;
        d = floor(linspace(1,ngauss,nslices+1));
        
        % do this in parallel !!!
        warning off
        delete(gcp('nocreate'))
        my_pp=parpool(ncores);
        parfor s = 1:nslices
            st = d(s);
            ed = d(s+1)-1;
            for gz = st:1:ed
                
                %Make design matrix for this predicted timecourse (no constant term!)
                %Fit the GLM using the matrix form
                %b = pinv(dm)*pixmat; %pinv = inv(dm'*dm)*dm'
                dm = [PRF_pred(:,gz),ones(size(PRF_pred,1),1)];
                [b,~,mss] = lscov(dm,MUAm); %get best fitting betas
                
                if ~Do.AllowNegBetas
                    %Set negative betas to Inf
                    mss(b(1,:)<0) = Inf; 
                end
                
                %Compare sse to current best, note down impriovements in teh minimum
                slice(s).ind(mss<slice(s).currentbest) = gz;
                slice(s).bestbeta(:,mss<slice(s).currentbest) = b(:,mss<slice(s).currentbest);
                slice(s).currentbest(mss<slice(s).currentbest) = mss(mss<slice(s).currentbest); %Update currentbest
            end
        end
        warning on
        
        %% Combine slices
        best = zeros(nslices,nchans);
        allbeta = zeros(2,nchans,nslices);
        allind = zeros(nslices,nchans);
        for s = 1:nslices
            best(s,:) = slice(s).currentbest;
            allbeta(:,:,s) = slice(s).bestbeta;
            allind(s,:) = slice(s).ind;
        end
        [bestval,bestix] = min(best);
        clear bestbeta
        clear ind
        for i = 1:nchans
            ind(i) = allind(bestix(i),i);
            bestbeta(:,i) = allbeta(:,i,bestix(i));
        end
        
        %READ in values
        rmx = gdets(ind,1);
        rmy = gdets(ind,2);
        rsd = gdets(ind,3);
        
        %Convert sd back into FWHM (keep both!!!)
        fwhm = 2*sqrt(log(2)*2);
        rsd_fwhm = rsd.*fwhm;
        
        %% Calculate corrcoef
        %turn into matrix computation
        rcrf = NaN(1,nchans);
        for N = 1:nchans
            %Get this channels
            pixel = MUAm(:,N);
            dm = [PRF_pred(:,ind(N)),ones(size(PRF_pred,1),1)]; %No offset allowed!
            %Fit the GLM using the matrix form
            b = pinv(dm)*pixel; %pinv = inv(dm'*dm)*dm'
            %Gives a best-fitting beta for each pixel.
            %Make the predicted response using the GLM equation y = x*beta
            y = dm*b;
            c = corrcoef(y,MUAm(:,N));
            rcrf(N) = c(1,2);
            
            bb = bestbeta(:,N);
            predTC(:,N) = PRF_pred(:,ind(N)).*bb(1)+bb(2);
        end
        
        %% Collect results
        PRFMUA.MUAm=MUAm;
        PRFMUA.rmx=rmx;    PRFMUA.rmy=rmy;
        PRFMUA.rsd=rsd;    PRFMUA.rcrf=rcrf;
        PRFMUA.rsd_fwhm=rsd_fwhm;
        PRFMUA.rsd_fwhm=rsd_fwhm;
        PRFMUA.beta=bestbeta;
        
        % collect all instances in 1 structure/file
        PRF_mua(a).MUAm=MUAm;
        PRF_mua(a).rmx=rmx;    PRF_mua(a).rmy=rmy;
        PRF_mua(a).rsd=rsd;    PRF_mua(a).rcrf=rcrf;
        PRF_mua(a).rsd_fwhm=rsd_fwhm;
        PRF_mua(a).pred=predTC;
        PRF_mua(a).beta=bestbeta;
        
        % save instance-specific
        fprintf('Saving results for instance %d\n', a);
        save(fullfile(save_fld,subj,sess,'MUA',...
            [subj '_' sess '_inst_' num2str(a) '_PRFMUA']),'PRFMUA','C');
        
        if 0
            %Example fit
            [i,j] = max(rcrf);
            %Get response from this pixel
            bestpix = MUAm(:,j);
            bb = bestbeta(:,j);
            %Get the predicted time-course for the best-fitting Gaussian
            %REMEMBER ind is the index of the best-fitting timecourse for each pixel
            predTC = PRF_pred(:,ind(j)).*bb(1)+bb(2);
            figure,bar(1:size(PRF_pred,1),bestpix),hold on,...
                scatter(1:size(PRF_pred,1),predTC,'r','filled')
            for sw=0:30:240
                plot([sw sw],[-0.5 3],'g');
            end
        end
        
        %     figure
        %     for j = 1:nchans
        %         bestpix = MUAm(:,j);
        %         bb = bestbeta(:,j);
        %         %Get the predicted time-course for the best-fitting Gaussian
        %         %REMEMBER ind is the index of the best-fitting timecourse for each pixel
        %         predTC = PRF_pred(:,ind(j)).*bb(1)+bb(2);
        %         bar(1:size(PRF_pred,1),bestpix),hold on,...
        %             scatter(1:size(PRF_pred,1),predTC,'r','filled')
        %         ylim([bb(2)-bb(1).*1.5 bb(2)+bb(1).*1.5])
        %         drawnow
        %         pause(0.1)
        %         hold off
        %     end
    end
    %% save
    save(fullfile(save_fld,subj,sess,'MUA',...
        [subj '_' sess '_allinst_PRFMUA']),'PRF_mua','C');
end

%% Load LFP responses & find best fitting pRFs ============================
if Do.FitLFP
    fprintf('Loading and processing LFP pRF results...\n');
    for a=1:8
        %% start
        fprintf('Array %d\n', a);
        load(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_array_' num2str(a) '_mLFP']));
        
        for fb=1:length(mLFP(a).freq)
            fprintf(['Freq' num2str(fb) ' ']);
            LFPm(fb).resp=[];
            for c=1:length(mLFP)
                LFPm(fb).resp=[LFPm(fb).resp ...
                    mLFP(c).freq(fb).bar'-mLFP(c).freq(fb).BL];
            end
            
            nchans = size(mLFP,2);
            for s = 1:12
                slice(s).ind = ones(nchans,1);
                slice(s).bestbeta = NaN(2,nchans);
                slice(s).currentbest = inf(1,nchans);
            end
            
            %% Divide the timecourse up into slices
            ngauss = size(PRF_pred,2);
            nslices = 12;
            d = floor(linspace(1,ngauss,nslices+1));
            
            % potentially do this in parallel !!!
            warning off
            delete(gcp('nocreate'))
            my_pp=parpool(ncores);
            parfor s = 1:nslices
                st = d(s);
                ed = d(s+1)-1;
                for gz = st:1:ed
                    
                    %Make design matrix for this predicted timecourse (no constant term!)
                    %Fit the GLM using the matrix form
                    %b = pinv(dm)*pixmat; %pinv = inv(dm'*dm)*dm'
                    dm = [PRF_pred(:,gz),ones(size(PRF_pred,1),1)];
                    [b,~,mss] = lscov(dm,LFPm(fb).resp); %get best fitting betas
 
                    if ~Do.AllowNegBetas
                        %Set negative betas to Inf
                        mss(b(1,:)<0) = Inf; 
                    end
                    
                    %Compare sse to current best, note down impriovements in teh minimum
                    slice(s).ind(mss<slice(s).currentbest) = gz;
                    slice(s).bestbeta(:,mss<slice(s).currentbest) = b(:,mss<slice(s).currentbest);
                    slice(s).currentbest(mss<slice(s).currentbest) = mss(mss<slice(s).currentbest); %Update currentbest
                end
            end
            warning on
            
            %% Combine slices
            best = zeros(nslices,nchans);
            allbeta = zeros(2,nchans,nslices);
            allind = zeros(nslices,nchans);
            for s = 1:nslices
                best(s,:) = slice(s).currentbest;
                allbeta(:,:,s) = slice(s).bestbeta;
                allind(s,:) = slice(s).ind;
            end
            [bestval,bestix] = min(best);
            clear bestbeta
            clear ind
            for i = 1:nchans
                ind(i) = allind(bestix(i),i);
                bestbeta(:,i) = allbeta(:,i,bestix(i));
            end
            
            %READ in values
            rmx = gdets(ind,1);
            rmy = gdets(ind,2);
            rsd = gdets(ind,3);
            
            %Convert sd back into FWHM
            fwhm = 2*sqrt(log(2)*2);
            rsd_fwhm = rsd.*fwhm;
            
            %% Calculate corrcoef
            %turn into matrix computation
            rcrf = NaN(1,nchans);
            for N = 1:nchans
                %Get this channels
                pixel = LFPm(fb).resp(:,N);
                dm = [PRF_pred(:,ind(N)),ones(size(PRF_pred,1),1)]; %No offset allowed!
                %Fit the GLM using the matrix form
                b = pinv(dm)*pixel; %pinv = inv(dm'*dm)*dm'
                %Gives a best-fitting beta for each pixel.
                %Make the predicted response using the GLM equation y = x*beta
                y = dm*b;
                c = corrcoef(y,LFPm(fb).resp(:,N));
                rcrf(N) = c(1,2);
                
                bb = bestbeta(:,N);
                predTC(:,N) = PRF_pred(:,ind(N)).*bb(1)+bb(2);
            end
            
            PRFLFP(fb).LFPm=LFPm(fb).resp;
            PRFLFP(fb).rmx=rmx;    PRFLFP(fb).rmy=rmy;
            PRFLFP(fb).rsd=rsd;    PRFLFP(fb).rcrf=rcrf;
            PRFLFP(fb).rsd_fwhm=rsd_fwhm;
            PRFLFP(fb).beta=bestbeta;
            PRFLFP(fb).pred=predTC;
            
            PRF_lfp(a,fb).LFPm=LFPm(fb).resp;
            PRF_lfp(a,fb).rmx=rmx;    PRF_lfp(a,fb).rmy=rmy;
            PRF_lfp(a,fb).rsd=rsd;    PRF_lfp(a,fb).rcrf=rcrf;
            PRF_lfp(a,fb).rsd_fwhm=rsd_fwhm;
            PRF_lfp(a,fb).beta=bestbeta;
            PRF_lfp(a,fb).pred=predTC;

        end
        
        fprintf('\nSaving results for instance %d\n', a);
        save(fullfile(save_fld,subj,sess,'LFP',...
            [subj '_' sess '_inst_' num2str(a) '_PRFLFP']),'PRFLFP','C');

    end
    save(fullfile(save_fld,subj,sess,'LFP',...
        [subj '_' sess '_allinst_PRFLFP']),'PRF_lfp','C');
end