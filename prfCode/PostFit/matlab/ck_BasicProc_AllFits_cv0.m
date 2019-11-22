function ck_BasicProc_AllFits_cv0

% this script will do some basic processing like calculate 
% - means from two-way crossvalidated results
% - X,Y coordinates from polar coordinates
% - highest R2 for crossval
% - distance between crossval runs in location and size

%%
% clear all;
clc;
CVMODE='cv0';

%% LOAD DATA ==============================================================
fitres_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/';

fprintf('==============================\n')
fprintf('Loading data...');
load(fullfile(fitres_path,'MultiModal',['AllFits_' CVMODE]),'R_MRI','R_EPHYS');
fprintf('DONE\n')
fprintf('==============================\n')

output_path = [...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/'...
    'Projects/pRF/FitResults/MultiModal/' CVMODE];
[~,~] = mkdir(output_path);

%% MRI --------------------------------------------------------------------
fprintf('MRI data ---------------------\n');
for r = 1:length(R_MRI) % animals
    for m = 1:length(R_MRI(r).model) % model fits
        fprintf(['Processing MODEL ' R_MRI(r).model(m).prfmodel '\n']);
        
        % convert pixels to degrees
        R_MRI(r).model(m).rfs = R_MRI(r).model(m).rfs./10; %#ok<*SAGROW>
        R_MRI(r).model(m).ecc = R_MRI(r).model(m).ecc./10;
        R_MRI(r).model(m).fwhm = R_MRI(r).model(m).fwhm./10;
                  
        % XY
        R_MRI(r).model(m).X = R_MRI(r).model(m).ecc.*cosd(R_MRI(r).model(m).ang);
        R_MRI(r).model(m).Y = R_MRI(r).model(m).ecc.*sind(R_MRI(r).model(m).ang);
        
        % angle
        rawang = atand(R_MRI(r).model(m).Y./R_MRI(r).model(m).X);
        rawang(R_MRI(r).model(m).X<0 & R_MRI(r).model(m).Y<0) = ...
            rawang(R_MRI(r).model(m).X<0 & R_MRI(r).model(m).Y<0)-180;
        rawang(R_MRI(r).model(m).X<0 & R_MRI(r).model(m).Y>=0) = ...
            rawang(R_MRI(r).model(m).X<0 & R_MRI(r).model(m).Y>=0)+180;
        R_MRI(r).model(m).ang = rawang;
        
    end
    % ROIs ----
    R_MRI(r).Brainmask = R_MRI(r).BRAIN>0;
    
    roi_idx = [];
    for rr = 1:length(R_MRI(r).ROI)
        roi_names{rr} = R_MRI(r).ROI(rr).label;
        roi_idx = [roi_idx R_MRI(r).ROI(rr).idx];
    end
    R_MRI(r).roi_names = roi_names;
    R_MRI(r).roi_idx = roi_idx;
end

%% EPHYS ------------------------------------------------------------------
fprintf('Ephys data ---------------------\n');
for r = 1:length(R_EPHYS) % animals
    for m = 1:length(R_EPHYS(r).model) % model fits
        fprintf(['Processing MODEL ' R_EPHYS(r).model(m).prfmodel '\n']);
        if ~strcmp(R_EPHYS(r).model(m).prfmodel,'classicRF')
            
            %% MUA --------------------------------------------------------
            fprintf('MUA ');
            for i = 1:8 % all instances
                % make name consistent with MRI
                R_EPHYS(r).model(m).MUA(i).rfs = R_EPHYS(r).model(m).MUA(i).rfsize;
                
                % convert pixels to degrees
                R_EPHYS(r).model(m).MUA(i).ecc = ...
                    R_EPHYS(r).model(m).MUA(i).ecc .* R_EPHYS(r).model(m).MUA(i).Pix2Deg;
                R_EPHYS(r).model(m).MUA(i).rfs = ...
                    R_EPHYS(r).model(m).MUA(i).rfs .* R_EPHYS(r).model(m).MUA(i).Pix2Deg;
                
                % calculate fwhm
                R_EPHYS(r).model(m).MUA(i).fwhm = ...
                    R_EPHYS(r).model(m).MUA(i).rfs.*(2*sqrt(2*log(2)));
                
                % XY
                R_EPHYS(r).model(m).MUA(i).X = ...
                    R_EPHYS(r).model(m).MUA(i).ecc.*cosd(R_EPHYS(r).model(m).MUA(i).ang);
                R_EPHYS(r).model(m).MUA(i).Y = ...
                    R_EPHYS(r).model(m).MUA(i).ecc.*sind(R_EPHYS(r).model(m).MUA(i).ang);
                
                % angle
                rawang = atand(R_EPHYS(r).model(m).MUA(i).Y./R_EPHYS(r).model(m).MUA(i).X);
                rawang(R_EPHYS(r).model(m).MUA(i).X<0 & R_EPHYS(r).model(m).MUA(i).Y<0) = ...
                    rawang(R_EPHYS(r).model(m).MUA(i).X<0 & R_EPHYS(r).model(m).MUA(i).Y<0)-180;
                rawang(R_EPHYS(r).model(m).MUA(i).X<0 & R_EPHYS(r).model(m).MUA(i).Y>=0) = ...
                    rawang(R_EPHYS(r).model(m).MUA(i).X<0 & R_EPHYS(r).model(m).MUA(i).Y>=0)+180;
                R_EPHYS(r).model(m).MUA(i).ang = rawang;
            end
            if isfield(R_EPHYS(r).model(m).MUA,'rfsize')
                R_EPHYS(r).model(m).MUA = rmfield(R_EPHYS(r).model(m).MUA,'rfsize');
            end
            
            %% LFP --------------------------------------------------------
            fprintf('LFP\n');
            for fb=1:5
                for i = 1:8 % all instances
                    % make name consistent with MRI
                    R_EPHYS(r).model(m).LFP(i,fb).rfs = R_EPHYS(r).model(m).LFP(i,fb).rfsize;
                    
                    % convert pixels to degrees
                    R_EPHYS(r).model(m).LFP(i,fb).ecc = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc .* R_EPHYS(r).model(m).LFP(i,fb).Pix2Deg;
                    R_EPHYS(r).model(m).LFP(i,fb).rfs = ...
                        R_EPHYS(r).model(m).LFP(i,fb).rfs .* R_EPHYS(r).model(m).LFP(i,fb).Pix2Deg;
                    
                    % calculate fwhm
                    R_EPHYS(r).model(m).LFP(i,fb).fwhm = ...
                        R_EPHYS(r).model(m).LFP(i,fb).rfs.*(2*sqrt(2*log(2)));
                
                    % XY
                    R_EPHYS(r).model(m).LFP(i,fb).X = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc.*cosd(R_EPHYS(r).model(m).LFP(i,fb).ang);
                    R_EPHYS(r).model(m).LFP(i,fb).Y = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc.*sind(R_EPHYS(r).model(m).LFP(i,fb).ang);
                    
                    % angle
                    rawang = atand(R_EPHYS(r).model(m).LFP(i,fb).Y./R_EPHYS(r).model(m).LFP(i,fb).X);
                    rawang(R_EPHYS(r).model(m).LFP(i,fb).X<0 & R_EPHYS(r).model(m).LFP(i,fb).Y<0) = ...
                        rawang(R_EPHYS(r).model(m).LFP(i,fb).X<0 & R_EPHYS(r).model(m).LFP(i,fb).Y<0)-180;
                    rawang(R_EPHYS(r).model(m).LFP(i,fb).X<0 & R_EPHYS(r).model(m).LFP(i,fb).Y>=0) = ...
                        rawang(R_EPHYS(r).model(m).LFP(i,fb).X<0 & R_EPHYS(r).model(m).LFP(i,fb).Y>=0)+180;
                    R_EPHYS(r).model(m).LFP(i,fb).ang = rawang;
                    
                end
            end
            if isfield(R_EPHYS(r).model(m).LFP,'rfsize')
                R_EPHYS(r).model(m).LFP = rmfield(R_EPHYS(r).model(m).LFP,'rfsize');
            end
            
        else % classic RF mapping
            fprintf('MUA\n');
            for i = 1:8 % instances
                for c = 1:length(R_EPHYS(r).model(m).MUA(i).RF) % channels
                    if c == 1 % pre-allocate
                        R_EPHYS(r).model(m).MUA(i).X = ...
                            nan(length(R_EPHYS(r).model(m).MUA(i).RF),1);
                        R_EPHYS(r).model(m).MUA(i).Y = ...
                            nan(length(R_EPHYS(r).model(m).MUA(i).RF),1);
                        R_EPHYS(r).model(m).MUA(i).rfs = ...
                            nan(length(R_EPHYS(r).model(m).MUA(i).RF),1);
                        R_EPHYS(r).model(m).MUA(i).ang = ...
                            nan(length(R_EPHYS(r).model(m).MUA(i).RF),1);
                        R_EPHYS(r).model(m).MUA(i).ecc = ...
                            nan(length(R_EPHYS(r).model(m).MUA(i).RF),1);
                    end
                    
                    % rfs / ang / ecc
                    R_EPHYS(r).model(m).MUA(i).rfs(c) = ...
                        R_EPHYS(r).model(m).MUA(i).RF{c}.szdeg;
                    R_EPHYS(r).model(m).MUA(i).ang(c) = ...
                        R_EPHYS(r).model(m).MUA(i).RF{c}.theta;
                    R_EPHYS(r).model(m).MUA(i).ecc(c) = ...
                        R_EPHYS(r).model(m).MUA(i).RF{c}.ecc;
                    
                    % XY
                    PixPerDeg = R_EPHYS(r).model(m).MUA(i).RF{c}.sz./...
                        R_EPHYS(r).model(m).MUA(i).RF{c}.szdeg;
                    R_EPHYS(r).model(m).MUA(i).X(c) = ...
                        R_EPHYS(r).model(m).MUA(i).RF{c}.centrex./PixPerDeg;
                    R_EPHYS(r).model(m).MUA(i).Y(c) = ...
                        R_EPHYS(r).model(m).MUA(i).RF{c}.centrey./PixPerDeg;

                    % fwhm
                    R_EPHYS(r).model(m).MUA(i).fwhm = ...
                        R_EPHYS(r).model(m).MUA(i).rfs.*(2*sqrt(2*log(2)));
                    
                    % SNR
                     R_EPHYS(r).model(m).MUA(i).SNR =  ...
                         R_EPHYS(r).model(m).MUA(i).SNR;
                end
            end
        end
    end 
    % ChannelMapping =====
    fprintf('Collecting the channel-map\n')
    CM = []; cm = R_EPHYS(r).ChanMap;
    for i = 1:size(cm.arrayNums,2) % instances
        for ch = 1:size(cm.arrayNums,1) % channels
            CM=[CM;...
                i ch ...
                cm.arrayNums(ch,i) cm.channelNums(ch,i) cm.areas(ch,i) ];
        end
    end
    R_EPHYS(r).cm = CM;
end

%% Create a bunch of tables with all results ==============================
% MRI ---
% Only brain voxels
% monkey mode model
fprintf('\n==============================\n')
fprintf('Creating combined tables...\n')
fprintf('==============================\n')

%% MRI --------------------------------------------------------------------
fprintf('MRI ==\n');
for r = 1:length(R_MRI) % animals
    % adjust ROI labels so they can be used as fields
    for ridx=1:length(R_MRI(r).ROI)
        R_MRI(r).ROI(ridx).label2=R_MRI(r).ROI(ridx).label;
        if ~isnan(str2double(R_MRI(r).ROI(ridx).label2(1)))
            R_MRI(r).ROI(ridx).label2 = ['a' R_MRI(r).ROI(ridx).label2];
        end
        if length(R_MRI(r).ROI(ridx).label2) >=11 && ...
                (strcmp(R_MRI(r).ROI(ridx).label2(1:11), 'Danny_LH_V1') || ... 
                strcmp(R_MRI(r).ROI(ridx).label2(1:10),'Eddy_LH_V1'))
            R_MRI(r).ROI(ridx).label2 = 'V1_electrodes';
        elseif length(R_MRI(r).ROI(ridx).label2) >=11 && ...
                (strcmp(R_MRI(r).ROI(ridx).label2(1:11), 'Danny_LH_V4') || ... 
                strcmp(R_MRI(r).ROI(ridx).label2(1:10), 'Eddy_LH_V4'))
            R_MRI(r).ROI(ridx).label2 = 'V4_electrodes';
        end
    end
    % brain voxels logical
    bm=R_MRI(r).Brainmask;
    
    %create tables
    for m = 1:length(R_MRI(r).model) % model fits  
        if r==1 && m==1
            % start the structures
            RTM.Monkey =[]; RTM.Mode =[]; RTM.Model =[];
            for ridx=1:length(R_MRI(r).ROI)
                RTM.(R_MRI(r).ROI(ridx).label2) = [];
            end
            RTM.R2 = []; RTM.rfs = []; RTM.fwhm = [];
            RTM.X = []; RTM.Y = [];
            RTM.ang = []; RTM.ecc = [];
            % not for all
            RTM.expt = []; RTM.sdratio = []; 
            RTM.gain = []; RTM.normamp = [];
        end
        nVox = sum(bm);
        
        % mean ====
        % labels
        RTM_Monkey = cell(nVox,1);
        for n=1:nVox; RTM_Monkey{n}=R_MRI(r).monkey;end
        RTM.Monkey = cat(1,RTM.Monkey,RTM_Monkey);
        %--
        RTM_Mode = cell(nVox,1);
        for n=1:nVox; RTM_Mode{n}=R_MRI(r).mode;end
        RTM.Mode = cat(1,RTM.Mode,RTM_Mode);
        %--
        RTM_Model = cell(nVox,1);
        for n=1:nVox; RTM_Model{n}=R_MRI(r).model(m).prfmodel;end 
        RTM.Model = cat(1,RTM.Model,RTM_Model);
        %--
        for ridx=1:length(R_MRI(r).ROI)
            RTM.(R_MRI(r).ROI(ridx).label2) = cat(1,...
                RTM.(R_MRI(r).ROI(ridx).label2),...
                R_MRI(r).ROI(ridx).idx(bm));
        end
        % values
        RTM.rfs = cat(1,RTM.rfs,R_MRI(r).model(m).rfs(bm)');
        RTM.fwhm = cat(1,RTM.fwhm,R_MRI(r).model(m).fwhm(bm)');
        RTM.X = cat(1,RTM.X,R_MRI(r).model(m).X(bm)');
        RTM.Y = cat(1,RTM.Y,R_MRI(r).model(m).Y(bm)');
        RTM.ang = cat(1,RTM.ang,R_MRI(r).model(m).ang(bm)');
        RTM.ecc = cat(1,RTM.ecc,R_MRI(r).model(m).ecc(bm)');
        
        % optional fields
        if isfield(R_MRI(r).model(m),'R2') && ~isempty(R_MRI(r).model(m).R2)
            RTM.R2 = cat(1,RTM.R2,R_MRI(r).model(m).R2(bm)');
        else
            RTM.R2 = cat(1,RTM.R2,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m),'gain') && ~isempty(R_MRI(r).model(m).gain)
            RTM.gain = cat(1,RTM.gain,R_MRI(r).model(m).gain(bm)');
        else
            RTM.gain = cat(1,RTM.gain,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m),'expt') && ~isempty(R_MRI(r).model(m).expt)
            RTM.expt = cat(1,RTM.expt,R_MRI(r).model(m).expt(bm)');
        else
            RTM.expt = cat(1,RTM.expt,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m),'sdratio') && ~isempty(R_MRI(r).model(m).sdratio)
            RTM.sdratio = cat(1,RTM.sdratio,R_MRI(r).model(m).sdratio(bm)');
        else
            RTM.sdratio = cat(1,RTM.sdratio,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m),'normamp') && ~isempty(R_MRI(r).model(m).normamp)
            RTM.normamp = cat(1,RTM.normamp,R_MRI(r).model(m).normamp(bm)');
        else
            RTM.normamp = cat(1,RTM.normamp,nan(nVox,1));
        end
    end
end


tMRI = struct2table(RTM);
MRI.RTE = RTM;

%% EPHYS MUA --------------------------------------------------------------
fprintf('EPHYS MUA ==\n');
for r = 1:length(R_EPHYS) % animals
    %create tables
    for m = 1:length(R_EPHYS(r).model) % model fits  
        if r==1 && m==1
            % start the structures
            % labels
            RTE.Monkey =[]; RTE.Mode = []; 
            RTE.Model = []; RTE.SigType = [];
            RTE.Array = []; RTE.Chan = []; RTE.Area = [];
            
            % results
            RTE.R2 = []; RTE.rfs = []; RTE.fwhm = [];
            RTE.X = []; RTE.Y = [];
            RTE.ang = []; RTE.ecc = [];
            
            % not for all
            RTE.expt =[]; RTE.sdratio = []; 
            RTE.gain = []; RTE.normamp = [];
            RTE.SNR = [];
        end
        
        nChan = size(R_EPHYS(1).model(1).MUA(1).ang,1);
        nInst = length(R_EPHYS(1).model(1).MUA);
        
        % MUA ---------
        % labels
        RTE_Monkey = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTE_Monkey{i} = R_EPHYS(r).monkey;
                i=i+1;
            end
        end
        RTE.Monkey = cat(1,RTE.Monkey,RTE_Monkey);
        %--
        RTE_Mode = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTE_Mode{i} = R_EPHYS(r).mode;
                i=i+1;
            end
        end
        RTE.Mode = cat(1,RTE.Mode,RTE_Mode);
        %--
        RTE_Model = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTE_Model{i} = R_EPHYS(r).model(m).prfmodel;
                i=i+1;
            end
        end
        RTE.Model = cat(1,RTE.Model,RTE_Model);
        %--
        RTE_SigType = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTE_SigType{i} = 'MUA';
                i=i+1;
            end
        end
        RTE.SigType = cat(1,RTE.SigType,RTE_SigType);
        %--
        if m==1
            RTE_Array = repmat(R_EPHYS(r).cm(:,3),length(R_EPHYS(r).model),1);
            RTE.Array = cat(1,RTE.Array,RTE_Array);
            RTE_Chan = repmat(R_EPHYS(r).cm(:,4),length(R_EPHYS(r).model),1);
            RTE.Chan = cat(1,RTE.Chan,RTE_Chan);
            RTE_Area = repmat(R_EPHYS(r).cm(:,5),length(R_EPHYS(r).model),1);
            RTE.Area = cat(1,RTE.Area,RTE_Area);
        end
        %-- 
        for i=1:nInst
            if m==5 % classic RF
                % standard fields
                RTE.rfs = cat(1,RTE.rfs,R_EPHYS(r).model(m).MUA(i).rfs);
                RTE.fwhm = cat(1,RTE.fwhm,R_EPHYS(r).model(m).MUA(i).fwhm);
                RTE.X = cat(1,RTE.X,R_EPHYS(r).model(m).MUA(i).X);
                RTE.Y = cat(1,RTE.Y,R_EPHYS(r).model(m).MUA(i).Y);
                RTE.ang = cat(1,RTE.ang,R_EPHYS(r).model(m).MUA(i).ang);
                RTE.ecc = cat(1,RTE.ecc,R_EPHYS(r).model(m).MUA(i).ecc);
                % optional fields
                RTE.R2 = cat(1,RTE.R2,nan(nChan,1));
                RTE.gain = cat(1,RTE.gain,nan(nChan,1));
                RTE.expt = cat(1,RTE.expt,nan(nChan,1));
                RTE.sdratio = cat(1,RTE.sdratio,nan(nChan,1));
                RTE.normamp = cat(1,RTE.normamp,nan(nChan,1));
                RTE.SNR = cat(1,RTE.SNR,R_EPHYS(r).model(m).MUA(i).SNR);
            else
                % standard fields
                RTE.rfs = cat(1,RTE.rfs,R_EPHYS(r).model(m).MUA(i).rfs);
                RTE.fwhm = cat(1,RTE.fwhm,R_EPHYS(r).model(m).MUA(i).fwhm);
                RTE.X = cat(1,RTE.X,R_EPHYS(r).model(m).MUA(i).X);
                RTE.Y = cat(1,RTE.Y,R_EPHYS(r).model(m).MUA(i).Y);
                RTE.ang = cat(1,RTE.ang,R_EPHYS(r).model(m).MUA(i).ang);
                RTE.ecc = cat(1,RTE.ecc,R_EPHYS(r).model(m).MUA(i).ecc);
                % optional fields
                if isfield(R_EPHYS(r).model(m).MUA(i),'R2')
                    RTE.R2 = cat(1,RTE.R2,R_EPHYS(r).model(m).MUA(i).R2);
                else
                    RTE.R2 = cat(1,RTE.R2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i),'gain')
                    RTE.gain = cat(1,RTE.gain,R_EPHYS(r).model(m).MUA(i).gain);
                else
                    RTE.gain = cat(1,RTE.gain,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i),'expt')
                    RTE.expt = cat(1,RTE.expt,R_EPHYS(r).model(m).MUA(i).expt);
                else
                    RTE.expt = cat(1,RTE.expt,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i),'sdratio')
                    RTE.sdratio = cat(1,RTE.sdratio,R_EPHYS(r).model(m).MUA(i).sdratio);
                else
                    RTE.sdratio = cat(1,RTE.sdratio,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i),'normamp')
                    RTE.normamp = cat(1,RTE.normamp,R_EPHYS(r).model(m).MUA(i).normamp);
                else
                    RTE.normamp = cat(1,RTE.normamp,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i),'SNR')
                    RTE.SNR = cat(1,RTE.SNR,R_EPHYS(r).model(m).MUA(i).SNR);
                else
                    RTE.SNR = cat(1,RTE.SNR,nan(nChan,1));
                end
            end
        end
    end

    % remove the 5th modeltype from label fields where you don't want them
    nRows = nInst*nChan;
    sp=(r-1)*length(R_EPHYS(r).model)*nRows;

end
tMUA = struct2table(RTE);
MUA.RTE = RTE;

%% EPHYS LFP --------------------------------------------------------------
LFPlabels = {'Theta','Alpha','Beta','lGamma','hGamma'};
fprintf('EPHYS LFP ==\n');
for r = 1:length(R_EPHYS) % animals
    %create tables
    for m = 1:length(R_EPHYS(r).model)-1 % model fits (skip the last)
        if r==1 && m==1
            % start the structures
            % labels
            RTE.Monkey =[]; RTE.Mode = []; 
            RTE.Model = []; RTE.SigType = [];
            RTE.Array = []; RTE.Chan = []; RTE.Area = [];
            
            % results
            RTE.R2 = []; RTE.rfs = []; RTE.fwhm = [];
            RTE.X = []; RTE.Y = [];
            RTE.ang = []; RTE.ecc = [];
            
            % not for all
            RTE.expt =[]; RTE.sdratio = []; 
            RTE.gain = []; RTE.normamp = [];
            RTE.SNR = [];
        end
        
        nChan = size(R_EPHYS(1).model(1).LFP(1).ang,1);
        nInst = size(R_EPHYS(1).model(1).LFP,1);
        nFB = size(R_EPHYS(1).model(1).LFP,2);
        nMod = length(R_EPHYS(r).model)-1;
        
        % LFP ---------
        % labels
        RTE_Monkey = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                for f=1:nFB
                    RTE_Monkey{i} = R_EPHYS(r).monkey;
                    i=i+1;
                end
            end
        end
        RTE.Monkey = cat(1,RTE.Monkey,RTE_Monkey);
        %--
        RTE_Mode = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                for f=1:nFB
                    RTE_Mode{i} = R_EPHYS(r).mode;
                    i=i+1;
                end
            end
        end
        RTE.Mode = cat(1,RTE.Mode,RTE_Mode);
        %--
        RTE_Model = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                for f=1:nFB
                    RTE_Model{i} = R_EPHYS(r).model(m).prfmodel;
                    i=i+1;
                end
            end
        end
        RTE.Model = cat(1,RTE.Model,RTE_Model);
        %--
        RTE_SigType = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                for f=1:nFB
                    RTE_SigType{i} = LFPlabels{f};
                    i=i+1;
                end
            end
        end
        RTE.SigType = cat(1,RTE.SigType,RTE_SigType);
        %--
        if m==1
            RTE_Array=[];
            RTE_Chan=[];
            RTE_Area=[];
            idx=1;
            for ni = 1:nInst
                for n=1:nChan
                    for f=1:nFB
                        RTE_Array =[RTE_Array; R_EPHYS(r).cm(idx,3)];
                        RTE_Chan = [RTE_Chan; R_EPHYS(r).cm(idx,4)];   
                        RTE_Area = [RTE_Area; R_EPHYS(r).cm(idx,5)];
                    end
                    idx=idx+1;
                end
            end
            RTE.Array = [RTE.Array; repmat(RTE_Array,nMod,1)];
            RTE.Chan = [RTE.Chan; repmat(RTE_Chan,nMod,1)];
            RTE.Area = [RTE.Area; repmat(RTE_Area,nMod,1)];
        end
        
        %-- 
        for i=1:nInst
            for f=1:nFB
                % standard fields
                RTE.rfs = cat(1,RTE.rfs,R_EPHYS(r).model(m).LFP(i,f).rfs);
                RTE.fwhm = cat(1,RTE.fwhm,R_EPHYS(r).model(m).LFP(i,f).fwhm);
                RTE.X = cat(1,RTE.X,R_EPHYS(r).model(m).LFP(i,f).X);
                RTE.Y = cat(1,RTE.Y,R_EPHYS(r).model(m).LFP(i,f).Y);
                RTE.ang = cat(1,RTE.ang,R_EPHYS(r).model(m).LFP(i,f).ang);
                RTE.ecc = cat(1,RTE.ecc,R_EPHYS(r).model(m).LFP(i,f).ecc);
                % optional fields
                if isfield(R_EPHYS(r).model(m).LFP(i,f),'R2')
                    RTE.R2 = cat(1,RTE.R2,R_EPHYS(r).model(m).LFP(i,f).R2);
                else
                    RTE.R2 = cat(1,RTE.R2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f),'gain')
                    RTE.gain = cat(1,RTE.gain,R_EPHYS(r).model(m).LFP(i,f).gain);
                else
                    RTE.gain = cat(1,RTE.gain,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f),'expt')
                    RTE.expt = cat(1,RTE.expt,R_EPHYS(r).model(m).LFP(i,f).expt);
                else
                    RTE.expt = cat(1,RTE.expt,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f),'sdratio')
                    RTE.sdratio = cat(1,RTE.sdratio,R_EPHYS(r).model(m).LFP(i,f).sdratio);
                else
                    RTE.sdratio = cat(1,RTE.sdratio,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f),'normamp')
                    RTE.normamp = cat(1,RTE.normamp,R_EPHYS(r).model(m).LFP(i,f).normamp);
                else
                    RTE.normamp = cat(1,RTE.normamp,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f),'SNR')
                    RTE.SNR = cat(1,RTE.SNR,R_EPHYS(r).model(m).LFP(i,f).SNR);
                else
                    RTE.SNR = cat(1,RTE.SNR,nan(nChan,1));
                end
            end
        end
    end
end
    
tLFP = struct2table(RTE);
LFP.RTE = RTE;

%% SAVE the tables & structs ==============================================
fprintf('\n==============================\n')
fprintf('Saving the results\n')
fprintf('==============================\n')

fprintf('Tables...\n');

save(...
    fullfile(output_path,'Tables'),...
    'tMRI','tMUA','tLFP','-v7.3');

fprintf('Structures...\n');

save(...
    fullfile(output_path,'MRI_Struct'),'MRI','-v7.3');  
save(...
    fullfile(output_path,'MUA_Struct'),'MUA','-v7.3');  
save(...
    fullfile(output_path,'LFP_Struct'),'LFP','-v7.3');  

fprintf('\nALL DONE!\n');
