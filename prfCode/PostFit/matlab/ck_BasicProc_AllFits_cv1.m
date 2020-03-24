function ck_BasicProc_AllFits_cv1

% this script will do some basic processing 

%%
% clear all;
clc;
CVMODE='cv1';
dataset='ORG'; % AVG / NEW / ORG

%% LOAD DATA ==============================================================
fitres_path = ...
    '/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF/FitResults/';

fprintf('==============================\n')
fprintf('Loading data...');
load(fullfile(fitres_path,'MultiModal',dataset,...
    ['AllFits_' CVMODE]),'R_MRI','R_EPHYS','D99');
fprintf('DONE\n')
fprintf('==============================\n')

output_path = fullfile(...
    ['/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/'...
    'Projects/pRF/FitResults/MultiModal'], dataset, CVMODE);
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
        
        F={'rfs','gain','R2','fwhm','sdratio','normamp','expt'};
        
        % avg some values
        for f=1:length(F)
            if ~isempty(R_MRI(r).model(m).(F{f}))
                R_MRI(r).model(m).avg.(F{f}) = mean(R_MRI(r).model(m).(F{f}),1);
            end
        end
        
        % XY
        R_MRI(r).model(m).X = R_MRI(r).model(m).ecc.*cosd(R_MRI(r).model(m).ang);
        R_MRI(r).model(m).Y = R_MRI(r).model(m).ecc.*sind(R_MRI(r).model(m).ang);
        
        % avg XY
        R_MRI(r).model(m).avg.X=mean(R_MRI(r).model(m).X,1);
        R_MRI(r).model(m).avg.Y=mean(R_MRI(r).model(m).Y,1);
        R_MRI(r).model(m).avg.ecc=sqrt( ...
            R_MRI(r).model(m).avg.X.^2 + R_MRI(r).model(m).avg.Y.^2);
        
        % avg angle
        rawang = atand(R_MRI(r).model(m).avg.Y./R_MRI(r).model(m).avg.X);
        rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y<0) = ...
            rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y<0)-180;
        rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y>=0) = ...
            rawang(R_MRI(r).model(m).avg.X<0 & R_MRI(r).model(m).avg.Y>=0)+180;
        R_MRI(r).model(m).avg.ang = rawang;
        
        % dXY
        R_MRI(r).model(m).diff.XY = sqrt(...
            diff(R_MRI(r).model(m).X,1,1).^2 + ...
            diff(R_MRI(r).model(m).Y,1,1).^2);
        
        % dRFS
        R_MRI(r).model(m).diff.rfs = diff(R_MRI(r).model(m).rfs,1,1);
        
        % get maximum R2 version of some values
        [R_MRI(r).model(m).max.R2, R_MRI(r).model(m).max.R2_idx] = ...
            max(R_MRI(r).model(m).R2,[],1);
        R_MRI(r).model(m).max.R2_idx = ...
            [R_MRI(r).model(m).max.R2_idx==1; R_MRI(r).model(m).max.R2_idx==2];
        
        F={'rfs','gain','ecc','ang','fwhm','sdratio','normamp','expt'};
        for f=1:length(F)
            if ~isempty(R_MRI(r).model(m).(F{f}))
                R_MRI(r).model(m).max.(F{f}) = ...
                    R_MRI(r).model(m).(F{f})(R_MRI(r).model(m).max.R2_idx)';
            end
        end
        R_MRI(r).model(m).max.X = R_MRI(r).model(m).X(R_MRI(r).model(m).max.R2_idx)';
        R_MRI(r).model(m).max.Y = R_MRI(r).model(m).Y(R_MRI(r).model(m).max.R2_idx)';
    end
    % ROIs ----
    R_MRI(r).Brainmask = R_MRI(r).BRAIN>0;
    R_MRI(r).D99 = D99;
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
                              
                F={'rfs','gain','R2','fwhm','sdratio','normamp','expt'};
                
                for f=1:length(F)
                    if isfield(R_EPHYS(r).model(m).MUA(i), F{f})
                        R_EPHYS(r).model(m).MUA(i).avg.(F{f}) = mean(R_EPHYS(r).model(m).MUA(i).(F{f}),2);
                    end
                end
                
                % XY
                R_EPHYS(r).model(m).MUA(i).X = ...
                    R_EPHYS(r).model(m).MUA(i).ecc.*cosd(R_EPHYS(r).model(m).MUA(i).ang);
                R_EPHYS(r).model(m).MUA(i).Y = ...
                    R_EPHYS(r).model(m).MUA(i).ecc.*sind(R_EPHYS(r).model(m).MUA(i).ang);
                
                % avg XY
                R_EPHYS(r).model(m).MUA(i).avg.X=mean(R_EPHYS(r).model(m).MUA(i).X,2);
                R_EPHYS(r).model(m).MUA(i).avg.Y=mean(R_EPHYS(r).model(m).MUA(i).Y,2);
                R_EPHYS(r).model(m).MUA(i).avg.ecc=sqrt( ...
                    R_EPHYS(r).model(m).MUA(i).avg.X.^2 + R_EPHYS(r).model(m).MUA(i).avg.Y.^2);
                
                % avg angle
                rawang = atand(R_EPHYS(r).model(m).MUA(i).avg.Y./R_EPHYS(r).model(m).MUA(i).avg.X);
                rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y<0) = ...
                    rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y<0)-180;
                rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y>=0) = ...
                    rawang(R_EPHYS(r).model(m).MUA(i).avg.X<0 & R_EPHYS(r).model(m).MUA(i).avg.Y>=0)+180;
                R_EPHYS(r).model(m).MUA(i).avg.ang = rawang;
                
                % dXY
                R_EPHYS(r).model(m).MUA(i).diff.XY = sqrt(...
                    diff(R_EPHYS(r).model(m).MUA(i).X,1,2).^2 + ...
                    diff(R_EPHYS(r).model(m).MUA(i).Y,1,2).^2);
                % dRFS
                R_EPHYS(r).model(m).MUA(i).diff.rfs = diff(R_EPHYS(r).model(m).MUA(i).rfs,1,2);
                
                % get maximum R2 version of some values
                [R_EPHYS(r).model(m).MUA(i).max.R2, R_EPHYS(r).model(m).MUA(i).max.R2_idx] = ...
                    max(R_EPHYS(r).model(m).MUA(i).R2,[],2);                
                R_EPHYS(r).model(m).MUA(i).max.R2_idx = ...
                     [R_EPHYS(r).model(m).MUA(i).max.R2_idx==1,...
                      R_EPHYS(r).model(m).MUA(i).max.R2_idx==2]';
                
                F={'rfs','gain','ecc','ang','fwhm','sdratio','normamp','expt'};
                for f=1:length(F)
                    if isfield(R_EPHYS(r).model(m).MUA(i), F{f}) && ~isempty(R_EPHYS(r).model(m).MUA(i).(F{f}))
                       TEMP = R_EPHYS(r).model(m).MUA(i).(F{f})'; 
                       R_EPHYS(r).model(m).MUA(i).max.(F{f}) = ...
                            TEMP(R_EPHYS(r).model(m).MUA(i).max.R2_idx);
                    end
                end
                TEMPX = R_EPHYS(r).model(m).MUA(i).X';
                R_EPHYS(r).model(m).MUA(i).max.X = TEMPX(R_EPHYS(r).model(m).MUA(i).max.R2_idx);
                TEMPY = R_EPHYS(r).model(m).MUA(i).Y';
                R_EPHYS(r).model(m).MUA(i).max.Y = TEMPY(R_EPHYS(r).model(m).MUA(i).max.R2_idx);
                
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
                
                    F={'rfs','gain','R2','fwhm','sdratio','normamp','expt'};
                    
                    for f=1:length(F)
                        if isfield(R_EPHYS(r).model(m).LFP(i,fb), F{f})
                            R_EPHYS(r).model(m).LFP(i,fb).avg.(F{f}) = mean(R_EPHYS(r).model(m).LFP(i,fb).(F{f}),2);
                        end
                    end
                    
                    % XY
                    R_EPHYS(r).model(m).LFP(i,fb).X = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc.*cosd(R_EPHYS(r).model(m).LFP(i,fb).ang);
                    R_EPHYS(r).model(m).LFP(i,fb).Y = ...
                        R_EPHYS(r).model(m).LFP(i,fb).ecc.*sind(R_EPHYS(r).model(m).LFP(i,fb).ang);
                    
                    % avg XY
                    R_EPHYS(r).model(m).LFP(i,fb).avg.X=mean(R_EPHYS(r).model(m).LFP(i,fb).X,2);
                    R_EPHYS(r).model(m).LFP(i,fb).avg.Y=mean(R_EPHYS(r).model(m).LFP(i,fb).Y,2);
                    R_EPHYS(r).model(m).LFP(i,fb).avg.ecc=sqrt( ...
                        R_EPHYS(r).model(m).LFP(i,fb).avg.X.^2 + R_EPHYS(r).model(m).LFP(i,fb).avg.Y.^2);
                    
                    % avg angle
                    rawang = atand(R_EPHYS(r).model(m).LFP(i,fb).avg.Y./R_EPHYS(r).model(m).LFP(i,fb).avg.X);
                    rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y<0) = ...
                        rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y<0)-180;
                    rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y>=0) = ...
                        rawang(R_EPHYS(r).model(m).LFP(i,fb).avg.X<0 & R_EPHYS(r).model(m).LFP(i,fb).avg.Y>=0)+180;
                    R_EPHYS(r).model(m).LFP(i,fb).avg.ang = rawang;
                    
                    % dXY
                    R_EPHYS(r).model(m).LFP(i,fb).diff.XY = sqrt(...
                        diff(R_EPHYS(r).model(m).LFP(i,fb).X,1,2).^2 + ...
                        diff(R_EPHYS(r).model(m).LFP(i,fb).Y,1,2).^2);
                    % dRFS
                    R_EPHYS(r).model(m).LFP(i,fb).diff.rfs = diff(R_EPHYS(r).model(m).LFP(i,fb).rfs,1,2);
                    
                    % get maximum R2 version of some values
                    [R_EPHYS(r).model(m).LFP(i,fb).max.R2, R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx] = ...
                        max(R_EPHYS(r).model(m).LFP(i,fb).R2,[],2);
                    R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx = ...
                        [R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx==1,...
                        R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx==2]';
                    
                    F={'rfs','gain','ecc','ang','fwhm','sdratio','normamp','expt'};
                    for f=1:length(F)
                        if isfield(R_EPHYS(r).model(m).LFP(i,fb), F{f}) && ~isempty(R_EPHYS(r).model(m).LFP(i,fb).(F{f}))
                            TEMP = R_EPHYS(r).model(m).LFP(i,fb).(F{f})';
                            R_EPHYS(r).model(m).LFP(i,fb).max.(F{f}) = ...
                                TEMP(R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx);
                        end
                    end
                    TEMPX = R_EPHYS(r).model(m).LFP(i,fb).X';
                    R_EPHYS(r).model(m).LFP(i,fb).max.X = TEMPX(R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx);
                    TEMPY = R_EPHYS(r).model(m).LFP(i,fb).Y';
                    R_EPHYS(r).model(m).LFP(i,fb).max.Y = TEMPY(R_EPHYS(r).model(m).LFP(i,fb).max.R2_idx);
                    
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
% == avg / max ==
% monkey mode model
fprintf('\n==============================\n')
fprintf('Creating combined tables...\n')
fprintf('==============================\n')

%% MRI --------------------------------------------------------------------
fprintf('MRI ==\n');
for r = 1:length(R_MRI) % animals
    
    % brain voxels logical
    bm=R_MRI(r).Brainmask;
    
    %create tables
    for m = 1:length(R_MRI(r).model) % model fits  
        if r==1 && m==1
            % start the structures
            RTMm.Monkey =[]; RTMm.Mode =[]; RTMm.Model =[];
            RTMm.ROI = [];
            RTMm.R2 = []; RTMm.rfs = []; RTMm.fwhm = [];
            RTMm.X = []; RTMm.Y = [];
            RTMm.ang = []; RTMm.ecc = [];
            % not for all
            RTMm.expt = []; RTMm.sdratio = []; 
            RTMm.gain = []; RTMm.normamp = [];
                      
            RTM.Monkey =[]; RTM.Mode =[]; RTM.Model =[];
            RTM.ROI = [];

            RTM.R2_1 = []; RTM.rfs_1 = []; RTM.fwhm_1 = [];
            RTM.X_1 = []; RTM.Y_1 = [];
            RTM.ang_1 = []; RTM.ecc_1 = [];
            RTM.R2_2 = []; RTM.rfs_2 = []; RTM.fwhm_2 = [];
            RTM.X_2 = []; RTM.Y_2 = [];
            RTM.ang_2 = []; RTM.ecc_2 = [];
            RTM.dXY=[]; RTM.dRFS = [];
            RTM.gain_1=[]; RTM.gain_2=[];
            RTM.expt_1 =[]; RTM.expt_2 =[]; 
            RTM.sdratio_1 = []; RTM.sdratio_2 = [];
            RTM.normamp_1 = []; RTM.normamp_2 = [];
            
        end
        nVox = sum(bm);
        
        % mean ====
        % labels
        RTMm_Monkey = cell(nVox,1);
        RTMm_Mode = cell(nVox,1);
        RTMm_Model = cell(nVox,1);
        for n=1:nVox
            RTMm_Monkey{n}=R_MRI(r).monkey;
            RTMm_Mode{n}=R_MRI(r).mode;
            RTMm_Model{n}=R_MRI(r).model(m).prfmodel;
        end
        RTMm.Monkey = cat(1,RTMm.Monkey,RTMm_Monkey);
        RTMm.Mode = cat(1,RTMm.Mode,RTMm_Mode);
        RTMm.Model = cat(1,RTMm.Model,RTMm_Model);

        % values
        RTMm.ROI = cat(1,RTMm.ROI,R_MRI(r).ROI(bm));
        RTMm.rfs = cat(1,RTMm.rfs,R_MRI(r).model(m).avg.rfs(bm)');
        RTMm.fwhm = cat(1,RTMm.fwhm,R_MRI(r).model(m).avg.fwhm(bm)');
        RTMm.X = cat(1,RTMm.X,R_MRI(r).model(m).avg.X(bm)');
        RTMm.Y = cat(1,RTMm.Y,R_MRI(r).model(m).avg.Y(bm)');
        RTMm.ang = cat(1,RTMm.ang,R_MRI(r).model(m).avg.ang(bm)');
        RTMm.ecc = cat(1,RTMm.ecc,R_MRI(r).model(m).avg.ecc(bm)');
        
        % optional fields
        if isfield(R_MRI(r).model(m).avg,'R2')
            RTMm.R2 = cat(1,RTMm.R2,R_MRI(r).model(m).avg.R2(bm)');
        else
            RTMm.R2 = cat(1,RTMm.R2,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m).avg,'gain')
            RTMm.gain = cat(1,RTMm.gain,R_MRI(r).model(m).avg.gain(bm)');
        else
            RTMm.gain = cat(1,RTMm.gain,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m).avg,'expt')
            RTMm.expt = cat(1,RTMm.expt,R_MRI(r).model(m).avg.expt(bm)');
        else
            RTMm.expt = cat(1,RTMm.expt,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m).avg,'sdratio')
            RTMm.sdratio = cat(1,RTMm.sdratio,R_MRI(r).model(m).avg.sdratio(bm)');
        else
            RTMm.sdratio = cat(1,RTMm.sdratio,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model(m).avg,'normamp')
            RTMm.normamp = cat(1,RTMm.normamp,R_MRI(r).model(m).avg.normamp(bm)');
        else
            RTMm.normamp = cat(1,RTMm.normamp,nan(nVox,1));
        end
        
        
        % diff ====
        % labels
        RTM.Monkey = RTMm.Monkey;
        RTM.Mode = RTMm.Mode;
        RTM.Model = RTMm.Model;

        % values
        RTM.ROI = cat(1,RTM.ROI,R_MRI(r).ROI(bm));
        RTM.rfs_1 = cat(1,RTM.rfs_1,R_MRI(r).model(m).rfs(1,bm)');
        RTM.rfs_2 = cat(1,RTM.rfs_2,R_MRI(r).model(m).rfs(2,bm)');
        RTM.fwhm_1 = cat(1,RTM.fwhm_1,R_MRI(r).model(m).fwhm(1,bm)');
        RTM.fwhm_2 = cat(1,RTM.fwhm_2,R_MRI(r).model(m).fwhm(2,bm)');
        RTM.X_1 = cat(1,RTM.X_1,R_MRI(r).model(m).X(1,bm)');
        RTM.X_2 = cat(1,RTM.X_2,R_MRI(r).model(m).X(2,bm)');
        RTM.Y_1 = cat(1,RTM.Y_1,R_MRI(r).model(m).Y(1,bm)');
        RTM.Y_2 = cat(1,RTM.Y_2,R_MRI(r).model(m).Y(2,bm)');
        RTM.ang_1 = cat(1,RTM.ang_1,R_MRI(r).model(m).ang(1,bm)');
        RTM.ang_2 = cat(1,RTM.ang_2,R_MRI(r).model(m).ang(2,bm)');
        RTM.ecc_1 = cat(1,RTM.ecc_1,R_MRI(r).model(m).ecc(1,bm)');
        RTM.ecc_2 = cat(1,RTM.ecc_2,R_MRI(r).model(m).ecc(2,bm)');
        
        RTM.dXY = cat(1,RTM.dXY,R_MRI(r).model(m).diff.XY(bm)');
        RTM.dRFS = cat(1,RTM.dRFS,R_MRI(r).model(m).diff.rfs(bm)');       
        
        % optional fields
        if isfield(R_MRI(r).model,'R2') && ~isempty(R_MRI(r).model(m).R2)
            RTM.R2_1 = cat(1,RTM.R2_1,R_MRI(r).model(m).R2(1,bm)');
            RTM.R2_2 = cat(1,RTM.R2_2,R_MRI(r).model(m).R2(2,bm)');
        else
            RTM.R2_1 = cat(1,RTM.R2_1,nan(nVox,1));
            RTM.R2_2 = cat(1,RTM.R2_2,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model,'gain') && ~isempty(R_MRI(r).model(m).gain)
            RTM.gain_1 = cat(1,RTM.gain_1,R_MRI(r).model(m).gain(1,bm)');
            RTM.gain_2 = cat(1,RTM.gain_2,R_MRI(r).model(m).gain(2,bm)');
        else
            RTM.gain_1 = cat(1,RTM.gain_1,nan(nVox,1));
            RTM.gain_2 = cat(1,RTM.gain_2,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model,'expt') && ~isempty(R_MRI(r).model(m).expt)
            RTM.expt_1 = cat(1,RTM.expt_1,R_MRI(r).model(m).expt(1,bm)');
            RTM.expt_2 = cat(1,RTM.expt_2,R_MRI(r).model(m).expt(2,bm)');
        else
            RTM.expt_1 = cat(1,RTM.expt_1,nan(nVox,1));
            RTM.expt_2 = cat(1,RTM.expt_2,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model,'sdratio') && ~isempty(R_MRI(r).model(m).sdratio)
            RTM.sdratio_1 = cat(1,RTM.sdratio_1,R_MRI(r).model(m).sdratio(1,bm)');
            RTM.sdratio_2 = cat(1,RTM.sdratio_2,R_MRI(r).model(m).sdratio(2,bm)');
        else
            RTM.sdratio_1 = cat(1,RTM.sdratio_1,nan(nVox,1));
            RTM.sdratio_2 = cat(1,RTM.sdratio_2,nan(nVox,1));
        end
        
        if isfield(R_MRI(r).model,'normamp') && ~isempty(R_MRI(r).model(m).normamp)
            RTM.normamp_1 = cat(1,RTM.normamp_1,R_MRI(r).model(m).normamp(1,bm)');
            RTM.normamp_2 = cat(1,RTM.normamp_2,R_MRI(r).model(m).normamp(2,bm)');
        else
            RTM.normamp_1 = cat(1,RTM.normamp_1,nan(nVox,1));
            RTM.normamp_2 = cat(1,RTM.normamp_2,nan(nVox,1));
        end
    end
end

tMRI = struct2table(RTM);
tMRI_mean = struct2table(RTMm);

MRI.RTE = RTM;
MRI.RTEm = RTMm;

%% EPHYS MUA --------------------------------------------------------------
fprintf('EPHYS MUA ==\n');
for r = 1:length(R_EPHYS) % animals
    %create tables
    for m = 1:length(R_EPHYS(r).model) % model fits  
        if r==1 && m==1
            % start the structures
            % labels
            RTEm.Monkey =[]; RTEm.Mode = []; 
            RTEm.Model = []; RTEm.SigType = [];
            RTEm.Array = []; RTEm.Chan = []; RTEm.Area = [];
            
            % results
            RTEm.R2 = []; RTEm.rfs = []; RTEm.fwhm = [];
            RTEm.X = []; RTEm.Y = [];
            RTEm.ang = []; RTEm.ecc = [];
            
            % not for all
            RTEm.expt =[]; RTEm.sdratio = []; 
            RTEm.gain = []; RTEm.normamp = [];
            RTEm.SNR = [];
            
            %
            RTE.R2_1 = []; RTE.rfs_1 = []; RTE.fwhm_1 = [];
            RTE.X_1 = []; RTE.Y_1 = [];
            RTE.ang_1 = []; RTE.ecc_1 = [];
            
            RTE.R2_2 = []; RTE.rfs_2 = []; RTE.fwhm_2 = [];
            RTE.X_2 = []; RTE.Y_2 = [];
            RTE.ang_2 = []; RTE.ecc_2 = [];
            
            RTE.Monkey =[]; RTE.Mode = []; 
            RTE.Model = []; RTE.SigType = [];
            RTE.Array = []; RTE.Chan = []; RTE.Area = [];
            
            RTE.gain_1 = []; RTE.gain_2 = [];
            RTE.expt_1 =[]; RTE.expt_2 =[]; 
            RTE.sdratio_1 = []; RTE.sdratio_2 = [];
            RTE.normamp_1 = []; RTE.normamp_2 = [];

            RTE.dXY=[]; RTE.dRFS = [];
            
            RTE.SNR = []; 
        end
        
        nChan = size(R_EPHYS(1).model(1).MUA(1).ang,1);
        nInst = length(R_EPHYS(1).model(1).MUA);
        
        % MUA ---------
        % labels
        RTEm_Monkey = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTEm_Monkey{i} = R_EPHYS(r).monkey;
                i=i+1;
            end
        end
        RTEm.Monkey = cat(1,RTEm.Monkey,RTEm_Monkey);
        %--
        RTEm_Mode = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTEm_Mode{i} = R_EPHYS(r).mode;
                i=i+1;
            end
        end
        RTEm.Mode = cat(1,RTEm.Mode,RTEm_Mode);
        %--
        RTEm_Model = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTEm_Model{i} = R_EPHYS(r).model(m).prfmodel;
                i=i+1;
            end
        end
        RTEm.Model = cat(1,RTEm.Model,RTEm_Model);
        %--
        RTEm_SigType = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                RTEm_SigType{i} = 'MUA';
                i=i+1;
            end
        end
        RTEm.SigType = cat(1,RTEm.SigType,RTEm_SigType);
        %--
        if m==1
            RTEm_Array = repmat(R_EPHYS(r).cm(:,3),length(R_EPHYS(r).model),1);
            RTEm.Array = cat(1,RTEm.Array,RTEm_Array);
            RTEm_Chan = repmat(R_EPHYS(r).cm(:,4),length(R_EPHYS(r).model),1);
            RTEm.Chan = cat(1,RTEm.Chan,RTEm_Chan);
            RTEm_Area = repmat(R_EPHYS(r).cm(:,5),length(R_EPHYS(r).model),1);
            RTEm.Area = cat(1,RTEm.Area,RTEm_Area);
        end
        %-- 
        for i=1:nInst
            if strcmp(R_EPHYS(r).model(m).prfmodel,'classicRF')
            %if m==5 % classic RF  
                % standard fields
                RTEm.rfs = cat(1,RTEm.rfs,R_EPHYS(r).model(m).MUA(i).rfs);
                RTEm.fwhm = cat(1,RTEm.fwhm,R_EPHYS(r).model(m).MUA(i).fwhm);
                RTEm.X = cat(1,RTEm.X,R_EPHYS(r).model(m).MUA(i).X);
                RTEm.Y = cat(1,RTEm.Y,R_EPHYS(r).model(m).MUA(i).Y);
                RTEm.ang = cat(1,RTEm.ang,R_EPHYS(r).model(m).MUA(i).ang);
                RTEm.ecc = cat(1,RTEm.ecc,R_EPHYS(r).model(m).MUA(i).ecc);
                % optional fields
                RTEm.R2 = cat(1,RTEm.R2,nan(nChan,1));
                RTEm.gain = cat(1,RTEm.gain,nan(nChan,1));
                RTEm.expt = cat(1,RTEm.expt,nan(nChan,1));
                RTEm.sdratio = cat(1,RTEm.sdratio,nan(nChan,1));
                RTEm.normamp = cat(1,RTEm.normamp,nan(nChan,1));
                RTEm.SNR = cat(1,RTEm.SNR,R_EPHYS(r).model(m).MUA(i).SNR);
            else
                % standard fields
                RTEm.rfs = cat(1,RTEm.rfs,R_EPHYS(r).model(m).MUA(i).avg.rfs);
                RTEm.fwhm = cat(1,RTEm.fwhm,R_EPHYS(r).model(m).MUA(i).avg.fwhm);
                RTEm.X = cat(1,RTEm.X,R_EPHYS(r).model(m).MUA(i).avg.X);
                RTEm.Y = cat(1,RTEm.Y,R_EPHYS(r).model(m).MUA(i).avg.Y);
                RTEm.ang = cat(1,RTEm.ang,R_EPHYS(r).model(m).MUA(i).avg.ang);
                RTEm.ecc = cat(1,RTEm.ecc,R_EPHYS(r).model(m).MUA(i).avg.ecc);
                % optional fields
                if isfield(R_EPHYS(r).model(m).MUA(i).avg,'R2')
                    RTEm.R2 = cat(1,RTEm.R2,R_EPHYS(r).model(m).MUA(i).avg.R2);
                else
                    RTEm.R2 = cat(1,RTEm.R2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i).avg,'gain')
                    RTEm.gain = cat(1,RTEm.gain,R_EPHYS(r).model(m).MUA(i).avg.gain);
                else
                    RTEm.gain = cat(1,RTEm.gain,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i).avg,'expt')
                    RTEm.expt = cat(1,RTEm.expt,R_EPHYS(r).model(m).MUA(i).avg.expt);
                else
                    RTEm.expt = cat(1,RTEm.expt,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i).avg,'sdratio')
                    RTEm.sdratio = cat(1,RTEm.sdratio,R_EPHYS(r).model(m).MUA(i).avg.sdratio);
                else
                    RTEm.sdratio = cat(1,RTEm.sdratio,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i).avg,'normamp')
                    RTEm.normamp = cat(1,RTEm.normamp,R_EPHYS(r).model(m).MUA(i).avg.normamp);
                else
                    RTEm.normamp = cat(1,RTEm.normamp,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i),'SNR')
                    RTEm.SNR = cat(1,RTEm.SNR,R_EPHYS(r).model(m).MUA(i).SNR);
                else
                    RTEm.SNR = cat(1,RTEm.SNR,nan(nChan,1));
                end
            end
        end
        
        
        % diff ====    
        %--
        for i=1:nInst % exclude classic RF here (no point in adding)
            if ~strcmp(R_EPHYS(r).model(m).prfmodel,'classicRF')
                % standard fields
                RTE.rfs_1 = cat(1,RTE.rfs_1,R_EPHYS(r).model(m).MUA(i).rfs(:,1));
                RTE.fwhm_1 = cat(1,RTE.fwhm_1,R_EPHYS(r).model(m).MUA(i).fwhm(:,1));
                RTE.X_1 = cat(1,RTE.X_1,R_EPHYS(r).model(m).MUA(i).X(:,1));
                RTE.Y_1 = cat(1,RTE.Y_1,R_EPHYS(r).model(m).MUA(i).Y(:,1));
                RTE.ang_1 = cat(1,RTE.ang_1,R_EPHYS(r).model(m).MUA(i).ang(:,1));
                RTE.ecc_1 = cat(1,RTE.ecc_1,R_EPHYS(r).model(m).MUA(i).ecc(:,1));
                %--
                RTE.rfs_2 = cat(1,RTE.rfs_2,R_EPHYS(r).model(m).MUA(i).rfs(:,2));
                RTE.fwhm_2 = cat(1,RTE.fwhm_2,R_EPHYS(r).model(m).MUA(i).fwhm(:,2));
                RTE.X_2 = cat(1,RTE.X_2,R_EPHYS(r).model(m).MUA(i).X(:,2));
                RTE.Y_2 = cat(1,RTE.Y_2,R_EPHYS(r).model(m).MUA(i).Y(:,2));
                RTE.ang_2 = cat(1,RTE.ang_2,R_EPHYS(r).model(m).MUA(i).ang(:,2));
                RTE.ecc_2 = cat(1,RTE.ecc_2,R_EPHYS(r).model(m).MUA(i).ecc(:,2));
                %--
                RTE.dXY = cat(1,RTE.dXY,R_EPHYS(r).model(m).MUA(i).diff.XY);
                RTE.dRFS = cat(1,RTE.dRFS,R_EPHYS(r).model(m).MUA(i).diff.rfs);
                %--
                
                % optional fields
                if isfield(R_EPHYS(r).model(m).MUA,'R2')
                    RTE.R2_1 = cat(1,RTE.R2_1,R_EPHYS(r).model(m).MUA(i).R2(:,1));
                    RTE.R2_2 = cat(1,RTE.R2_2,R_EPHYS(r).model(m).MUA(i).R2(:,2));
                else
                    RTE.R2_1 = cat(1,RTE.R2_1,nan(nChan,1));
                    RTE.R2_2 = cat(1,RTE.R2_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA,'gain')
                    RTE.gain_1 = cat(1,RTE.gain_1,R_EPHYS(r).model(m).MUA(i).gain(:,1));
                    RTE.gain_2 = cat(1,RTE.gain_2,R_EPHYS(r).model(m).MUA(i).gain(:,2));
                else
                    RTE.gain_1 = cat(1,RTE.gain_1,nan(nChan,1));
                    RTE.gain_2 = cat(1,RTE.gain_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA,'expt')
                    RTE.expt_1 = cat(1,RTE.expt_1,R_EPHYS(r).model(m).MUA(i).expt(:,1));
                    RTE.expt_2 = cat(1,RTE.expt_2,R_EPHYS(r).model(m).MUA(i).expt(:,2));
                else
                    RTE.expt_1 = cat(1,RTE.expt_1,nan(nChan,1));
                    RTE.expt_2 = cat(1,RTE.expt_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA,'sdratio')
                    RTE.sdratio_1 = cat(1,RTE.sdratio_1,R_EPHYS(r).model(m).MUA(i).sdratio(:,1));
                    RTE.sdratio_2 = cat(1,RTE.sdratio_2,R_EPHYS(r).model(m).MUA(i).sdratio(:,2));
                else
                    RTE.sdratio_1 = cat(1,RTE.sdratio_1,nan(nChan,1));
                    RTE.sdratio_2 = cat(1,RTE.sdratio_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA,'normamp')
                    RTE.normamp_1 = cat(1,RTE.normamp_1,R_EPHYS(r).model(m).MUA(i).normamp(:,1));
                    RTE.normamp_2 = cat(1,RTE.normamp_2,R_EPHYS(r).model(m).MUA(i).normamp(:,2));
                else
                    RTE.normamp_1 = cat(1,RTE.normamp_1,nan(nChan,1));
                    RTE.normamp_2 = cat(1,RTE.normamp_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).MUA(i),'SNR')
                    RTE.SNR = cat(1,RTE.SNR,R_EPHYS(r).model(m).MUA(i).SNR);
                else
                    RTE.SNR = cat(1,RTE.SNR,nan(nChan,1));
                end
            end
        end
    end

    % remove the last modeltype from label fields where you don't want them
    nRows = nInst*nChan;
    sp=(r-1)*length(R_EPHYS(r).model)*nRows;
    
    RTE.Monkey = cat(1,RTE.Monkey,RTEm.Monkey(sp+1:end-nRows,:));
    RTE.Mode = cat(1,RTE.Mode,RTEm.Mode(sp+1:end-nRows,:));
    RTE.Model = cat(1,RTE.Model,RTEm.Model(sp+1:end-nRows,:));
    RTE.SigType = cat(1,RTE.SigType,RTEm.SigType(sp+1:end-nRows,:));
    RTE.Array = cat(1,RTE.Array,RTEm.Array(sp+1:end-nRows,:));
    RTE.Chan = cat(1,RTE.Chan,RTEm.Chan(sp+1:end-nRows,:));
    RTE.Area = cat(1,RTE.Area,RTEm.Area(sp+1:end-nRows,:));

end
tMUA = struct2table(RTE);
tMUA_mean = struct2table(RTEm);

MUA.RTE = RTE;
MUA.RTEm = RTEm;

%% EPHYS LFP --------------------------------------------------------------
LFPlabels = {'Theta','Alpha','Beta','lGamma','hGamma'};
fprintf('EPHYS LFP ==\n');
for r = 1:length(R_EPHYS) % animals
    %create tables
    for m = 1:length(R_EPHYS(r).model)-1 % model fits (skip the last)
        if r==1 && m==1
            % start the structures
            % labels
            RTEm.Monkey =[]; RTEm.Mode = []; 
            RTEm.Model = []; RTEm.SigType = [];
            RTEm.Array = []; RTEm.Chan = []; RTEm.Area = [];
            
            % results
            RTEm.R2 = []; RTEm.rfs = []; RTEm.fwhm = [];
            RTEm.X = []; RTEm.Y = [];
            RTEm.ang = []; RTEm.ecc = [];
            
            % not for all
            RTEm.expt =[]; RTEm.sdratio = []; 
            RTEm.gain = []; RTEm.normamp = [];
            RTEm.SNR = [];
            
            RTE.R2_1 = []; RTE.rfs_1 = []; RTE.fwhm_1 = [];
            RTE.X_1 = []; RTE.Y_1 = [];
            RTE.ang_1 = []; RTE.ecc_1 = [];
            
            RTE.R2_2 = []; RTE.rfs_2 = []; RTE.fwhm_2 = [];
            RTE.X_2 = []; RTE.Y_2 = [];
            RTE.ang_2 = []; RTE.ecc_2 = [];
                        
            RTE.Monkey =[]; RTE.Mode = []; 
            RTE.Model = []; RTE.SigType = [];
            RTE.Array = []; RTE.Chan = []; RTE.Area = [];
            
            RTE.gain_1 = []; RTE.gain_2 = [];
            RTE.expt_1 =[]; RTE.expt_2 =[]; 
            RTE.sdratio_1 = []; RTE.sdratio_2 = [];
            RTE.normamp_1 = []; RTE.normamp_2 = [];

            RTE.dXY=[]; RTE.dRFS = [];
            
            RTE.SNR = []; 
        end
        
        nChan = size(R_EPHYS(1).model(1).LFP(1).ang,1);
        nInst = size(R_EPHYS(1).model(1).LFP,1);
        nFB = size(R_EPHYS(1).model(1).LFP,2);
        nMod = length(R_EPHYS(r).model)-1;
        
        % LFP ---------
        % labels
        RTEm_Monkey = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                for f=1:nFB
                    RTEm_Monkey{i} = R_EPHYS(r).monkey;
                    i=i+1;
                end
            end
        end
        RTEm.Monkey = cat(1,RTEm.Monkey,RTEm_Monkey);
        %--
        RTEm_Mode = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                for f=1:nFB
                    RTEm_Mode{i} = R_EPHYS(r).mode;
                    i=i+1;
                end
            end
        end
        RTEm.Mode = cat(1,RTEm.Mode,RTEm_Mode);
        %--
        RTEm_Model = cell(nChan*nInst,1); i=1;
        for ni = 1:nInst
            for n=1:nChan
                for f=1:nFB
                    RTEm_Model{i} = R_EPHYS(r).model(m).prfmodel;
                    i=i+1;
                end
            end
        end
        RTEm.Model = cat(1,RTEm.Model,RTEm_Model);
        %--
        RTEm_SigType = cell(nChan*nInst,1); i=1;
        for f=1:nFB
            for ni = 1:nInst
                for n=1:nChan
                    RTEm_SigType{i} = LFPlabels{f};
                    i=i+1;
                end
            end
        end
        RTEm.SigType = cat(1,RTEm.SigType,RTEm_SigType);
        %--
        if m==1
            RTEm_Array=[];
            RTEm_Chan=[];
            RTEm_Area=[];
            for f=1:nFB
                RTEm_Array=[RTEm_Array; R_EPHYS(r).cm(:,3)];
                RTEm_Chan=[RTEm_Chan; R_EPHYS(r).cm(:,4)];
                RTEm_Area=[RTEm_Area; R_EPHYS(r).cm(:,5)];
            end
            RTEm.Array = [RTEm.Array; repmat(RTEm_Array,nMod,1)];
            RTEm.Chan = [RTEm.Chan; repmat(RTEm_Chan,nMod,1)];
            RTEm.Area = [RTEm.Area; repmat(RTEm_Area,nMod,1)];
        end
        
        %-- 
        for f=1:nFB
            for i=1:nInst
                % standard fields
                RTEm.rfs = cat(1,RTEm.rfs,R_EPHYS(r).model(m).LFP(i,f).avg.rfs);
                RTEm.fwhm = cat(1,RTEm.fwhm,R_EPHYS(r).model(m).LFP(i,f).avg.fwhm);
                RTEm.X = cat(1,RTEm.X,R_EPHYS(r).model(m).LFP(i,f).avg.X);
                RTEm.Y = cat(1,RTEm.Y,R_EPHYS(r).model(m).LFP(i,f).avg.Y);
                RTEm.ang = cat(1,RTEm.ang,R_EPHYS(r).model(m).LFP(i,f).avg.ang);
                RTEm.ecc = cat(1,RTEm.ecc,R_EPHYS(r).model(m).LFP(i,f).avg.ecc);
                % optional fields
                if isfield(R_EPHYS(r).model(m).LFP(i,f).avg,'R2')
                    RTEm.R2 = cat(1,RTEm.R2,R_EPHYS(r).model(m).LFP(i,f).avg.R2);
                else
                    RTEm.R2 = cat(1,RTEm.R2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f).avg,'gain')
                    RTEm.gain = cat(1,RTEm.gain,R_EPHYS(r).model(m).LFP(i,f).avg.gain);
                else
                    RTEm.gain = cat(1,RTEm.gain,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f).avg,'expt')
                    RTEm.expt = cat(1,RTEm.expt,R_EPHYS(r).model(m).LFP(i,f).avg.expt);
                else
                    RTEm.expt = cat(1,RTEm.expt,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f).avg,'sdratio')
                    RTEm.sdratio = cat(1,RTEm.sdratio,R_EPHYS(r).model(m).LFP(i,f).avg.sdratio);
                else
                    RTEm.sdratio = cat(1,RTEm.sdratio,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f).avg,'normamp')
                    RTEm.normamp = cat(1,RTEm.normamp,R_EPHYS(r).model(m).LFP(i,f).avg.normamp);
                else
                    RTEm.normamp = cat(1,RTEm.normamp,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP(i,f),'SNR')
                    RTEm.SNR = cat(1,RTEm.SNR,R_EPHYS(r).model(m).LFP(i,f).SNR);
                else
                    RTEm.SNR = cat(1,RTEm.SNR,nan(nChan,1));
                end
            end
        end
        
        
        % diff ====    
        %--
        for f=1:nFB
            for i=1:nInst 
                % standard fields
                RTE.rfs_1 = cat(1,RTE.rfs_1,R_EPHYS(r).model(m).LFP(i,f).rfs(:,1));
                RTE.fwhm_1 = cat(1,RTE.fwhm_1,R_EPHYS(r).model(m).LFP(i,f).fwhm(:,1));
                RTE.X_1 = cat(1,RTE.X_1,R_EPHYS(r).model(m).LFP(i,f).X(:,1));
                RTE.Y_1 = cat(1,RTE.Y_1,R_EPHYS(r).model(m).LFP(i,f).Y(:,1));
                RTE.ang_1 = cat(1,RTE.ang_1,R_EPHYS(r).model(m).LFP(i,f).ang(:,1));
                RTE.ecc_1 = cat(1,RTE.ecc_1,R_EPHYS(r).model(m).LFP(i,f).ecc(:,1));
                %--
                RTE.rfs_2 = cat(1,RTE.rfs_2,R_EPHYS(r).model(m).LFP(i,f).rfs(:,2));
                RTE.fwhm_2 = cat(1,RTE.fwhm_2,R_EPHYS(r).model(m).LFP(i,f).fwhm(:,2));
                RTE.X_2 = cat(1,RTE.X_2,R_EPHYS(r).model(m).LFP(i,f).X(:,2));
                RTE.Y_2 = cat(1,RTE.Y_2,R_EPHYS(r).model(m).LFP(i,f).Y(:,2));
                RTE.ang_2 = cat(1,RTE.ang_2,R_EPHYS(r).model(m).LFP(i,f).ang(:,2));
                RTE.ecc_2 = cat(1,RTE.ecc_2,R_EPHYS(r).model(m).LFP(i,f).ecc(:,2));
                %--
                RTE.dXY = cat(1,RTE.dXY,R_EPHYS(r).model(m).LFP(i,f).diff.XY);
                RTE.dRFS = cat(1,RTE.dRFS,R_EPHYS(r).model(m).LFP(i,f).diff.rfs);
                %--
                
                % optional fields
                if isfield(R_EPHYS(r).model(m).LFP,'R2')
                    RTE.R2_1 = cat(1,RTE.R2_1,R_EPHYS(r).model(m).LFP(i,f).R2(:,1));
                    RTE.R2_2 = cat(1,RTE.R2_2,R_EPHYS(r).model(m).LFP(i,f).R2(:,2));
                else
                    RTE.R2_1 = cat(1,RTE.R2_1,nan(nChan,1));
                    RTE.R2_2 = cat(1,RTE.R2_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP,'gain')
                    RTE.gain_1 = cat(1,RTE.gain_1,R_EPHYS(r).model(m).LFP(i,f).gain(:,1));
                    RTE.gain_2 = cat(1,RTE.gain_2,R_EPHYS(r).model(m).LFP(i,f).gain(:,2));
                else
                    RTE.gain_1 = cat(1,RTE.gain_1,nan(nChan,1));
                    RTE.gain_2 = cat(1,RTE.gain_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP,'expt')
                    RTE.expt_1 = cat(1,RTE.expt_1,R_EPHYS(r).model(m).LFP(i,f).expt(:,1));
                    RTE.expt_2 = cat(1,RTE.expt_2,R_EPHYS(r).model(m).LFP(i,f).expt(:,2));
                else
                    RTE.expt_1 = cat(1,RTE.expt_1,nan(nChan,1));
                    RTE.expt_2 = cat(1,RTE.expt_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP,'sdratio')
                    RTE.sdratio_1 = cat(1,RTE.sdratio_1,R_EPHYS(r).model(m).LFP(i,f).sdratio(:,1));
                    RTE.sdratio_2 = cat(1,RTE.sdratio_2,R_EPHYS(r).model(m).LFP(i,f).sdratio(:,2));
                else
                    RTE.sdratio_1 = cat(1,RTE.sdratio_1,nan(nChan,1));
                    RTE.sdratio_2 = cat(1,RTE.sdratio_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP,'normamp')
                    RTE.normamp_1 = cat(1,RTE.normamp_1,R_EPHYS(r).model(m).LFP(i,f).normamp(:,1));
                    RTE.normamp_2 = cat(1,RTE.normamp_2,R_EPHYS(r).model(m).LFP(i,f).normamp(:,2));
                else
                    RTE.normamp_1 = cat(1,RTE.normamp_1,nan(nChan,1));
                    RTE.normamp_2 = cat(1,RTE.normamp_2,nan(nChan,1));
                end
                
                if isfield(R_EPHYS(r).model(m).LFP,'SNR')
                    RTE.SNR = cat(1,RTE.SNR,R_EPHYS(r).model(m).LFP(i,f).SNR);
                else
                    RTE.SNR = cat(1,RTE.SNR,nan(nChan,1));
                end
            end
        end
    end
end
RTE.Monkey = RTEm.Monkey;
RTE.Mode = RTEm.Mode;
RTE.Model = RTEm.Model;
RTE.SigType = RTEm.SigType;
RTE.Array = RTEm.Array;
RTE.Chan = RTEm.Chan;
RTE.Area = RTEm.Area;
    
tLFP = struct2table(RTE);
tLFP_mean = struct2table(RTEm);

LFP.RTE = RTE;
LFP.RTEm = RTEm;

%% SAVE the tables & structs ==============================================
fprintf('\n==============================\n')
fprintf('Saving the results\n')
fprintf('==============================\n')

fprintf('Tables...\n');
% save(...
%     fullfile(output_path,'CombiTables'),...
%     'tMRI','tMRI_mean',... 
%     'tMUA','tMUA_mean',... 
%     'tLFP','tLFP_mean','-v7.3');

save(...
    fullfile(output_path,'Tables_diff'),...
    'tMRI','tMUA','tLFP','-v7.3');
save(...
    fullfile(output_path,'Tables_mean'),...
    'tMRI_mean','tMUA_mean','tLFP_mean','-v7.3');

fprintf('Structures...\n');
% save(...
%     fullfile(output_path,'CombiStructs'),...
%     'MRI','MUA','LFP','-v7.3');    

save(...
    fullfile(output_path,'MRI_Struct'),'MRI','-v7.3');  
save(...
    fullfile(output_path,'MUA_Struct'),'MUA','-v7.3');  
save(...
    fullfile(output_path,'LFP_Struct'),'LFP','-v7.3');  

fprintf('\nALL DONE!\n');
