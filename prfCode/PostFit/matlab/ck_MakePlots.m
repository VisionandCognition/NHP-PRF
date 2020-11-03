% ck_MakePlots
% Takes the pre-processed fitting result tables and creates comparison plots

SaveFigs=false;
CloseFigs=false;

fprintf('Moving to root folder for this script\n')
% go to the location of this files
cd(fullfile('/Users','chris','Documents','MRI_ANALYSIS',...
    'NHP-PRF','prfCode','PostFit','matlab'));


%% ========================================================================
%  INITIATE ---------------------------------------------------------------
%  ========================================================================

%% Paths ==================================================================
fprintf('Setting up paths\n')
BaseFld = pwd;
DS='ORG';
ResFld = ...
    ['/Users/chris/Documents/MRI_ANALYSIS/NHP-PRF/'...
    'FitResults/MultiModal/' DS '/cv1'];
TT='Tables_mean';

% add colorbrewer
addpath(genpath('~/Documents/MATLAB/GENERAL_FUNCTIONS/BrewerMap'));
def_cmap = 'Spectral';
% add matplotlib
addpath(genpath('~/Documents/MATLAB/GENERAL_FUNCTIONS/matplotlib'));
% add boundedline
addpath(genpath('~/Documents/MATLAB/GENERAL_FUNCTIONS/boundedline'));
ResType = 'mean'; % max / mean
TT = ['Tables_' ResType];

% figure saving folder
pngfld = fullfile(pwd,'fig_png');
svgfld = fullfile(pwd,'fig_svg');
[~,~] = mkdir(pngfld); [~,~] = mkdir(svgfld);

%% Load ===================================================================
fprintf(['Loading results table. '...
    'Please be patient, this will take a while..\n']);
tic; load(fullfile(ResFld,TT)); t_load=toc;
fprintf(['Loading took ' num2str(t_load) ' s\n']);

%% Correct for result type ================================================
fprintf('Generalizing variable names...\n');
fprintf('MRI...'); 
eval(['tMRI = tMRI_' ResType '; clear tMRI_' ResType ';']);
fprintf('MUA...'); 
eval(['tMUA = tMUA_' ResType '; clear tMUA_' ResType ';']);
fprintf('LFP...'); 
eval(['tLFP = tLFP_' ResType '; clear tLFP_' ResType ';']);
fprintf('>> DONE\n');

%% Correct the eccentricity for Classic-RF from pix to deg ================
fprintf('Correcting eccentricity values for ClassicRF from pix to deg\n')
CRF_idx=strcmp(tMUA.Model,'classicRF');
tMUA.ecc(CRF_idx) = sqrt(tMUA.X(CRF_idx).^2 + tMUA.Y(CRF_idx).^2);

%% Correct the rf size for Classic-RF from diameter to radius =============
fprintf('Correcting RF size values for ClassicRF from diam to radius\n')
tMUA.rfs(CRF_idx) = tMUA.rfs(CRF_idx)./2;

%% Split MRI table by model ===============================================
fprintf('Splitting MRI table by model\n')
m = unique(tMRI.Model);
for mi = 1:length(m)
    M = tMRI(strcmp(tMRI.Model,m{mi}),:);
    T(mi).mod = M;
    T(mi).name = m{mi};
    modidx.(m{mi}) = mi;
end

%% Initiate some information ==============================================
fprintf('Generating ROI, model, and subject labels\n')
% proces ROI info
rois = {...
    'V1',   [34];...            % occipital / visual
    'V2',   [131];...           % occipital / visual
    'V3',   [60,93];...     	% occipital / visual
    'V3A'	[123];...			% occipital / visual
    'V4',   [20,39,75];...      % mid-visual
    'V6',   [73];...            % mid-visual
    'V6A',  [141,56];...
    'MT',   [95];...            % mid-visual
    'MST',  [99];...            % mid-visual
    'TEO',  [125];...           % temporal
    'TAa',  [152];...           % temporal
    'Tpt',  [97];...            % temporal (temporal/parietal)
    'TPO',  [159];...           % temporal
    'FST',  [53];...            % temporal
    'A1',	[80];...			% temporal
    'ML'	[25];...			% temporal	
    'AL',	[46];...			% temporal
    'PULV', [197,198];...       % subcortical
    'LGN',  [200,201];...       % subcortical
    'STR',  [175];...           % subcortical
    'MIP',  [89];...  
    'PIP',  [90];...
    'LIP',  [31,130];...        % parietal
    'VIP',  [30];...            % parietal
    '5',    [134];...         	% parietal
    '7',    [91,121];...        % parietal
    'SI',   [50,137];...        % Prim Somatosensory
    'SII',  [63];...            % Sec Somatosensory
    'F2',   [153];...           % premotor
    'F4',   [146];...           % premotor
    'F5',   [129];...           % premotor
    'F7',   [126];...           % premotor
    '8',    [32,51,57,148];...  % frontal FEF  8A: 51, 148; 8B: 32, 57 % combine these??
    'CINp', [6,17,27];...       % posterior cingulate
    'CINa', [45,98];...         % anterior cingulate
    'OFC',  [55,107];...        % frontal
    'INS',  [18,87,128];...     % frontal
    'DLPFC',[47,76,127];...     % frontal
    'VMPFC',[5];...         	% frontal
};

roi=rois(:,2);
roilabels=rois(:,1);

MRI_MODELS={...
    'linhrf_cv1_mhrf','linhrf_cv1_dhrf';...
    'linhrf_cv1_mhrf_neggain','linhrf_cv1_dhrf_neggain';...
    'csshrf_cv1_mhrf','csshrf_cv1_dhrf';...
    'doghrf_cv1_mhrf','doghrf_cv1_dhrf';...
    };
ephys_MOD={'linear_ephys_cv1','linear_ephys_cv1_neggain',...
    'css_ephys_cv1','dog_ephys_cv1'};
MMS={...
    'LIN_m','LIN_d';...
    'LIN-N_m','LIN-N_d';...
    'CSS_m','CSS_d';...
    'DOG_m','DOG_d';...
    };

SUBS = unique(tMRI.Monkey);

%% Inspect DoG size indication ============================================
fprintf('Handling DoG pRF sizes MRI \n')
R2th=20;
DoGtest = tMRI(...
    strcmp(tMRI.Model,'doghrf_cv1_mhrf') & ...
    tMRI.ROI == ck_GetROIidx({'V1'},rois) & ...
    tMRI.R2 > R2th ,:);

ftemp = figure;
subplot(2,2,1);
scatter(DoGtest.ecc,DoGtest.rfs);
xlabel('ecc');ylabel('rfs');
subplot(2,2,2);
scatter(DoGtest.rfs,DoGtest.gain)
xlabel('rfs');ylabel('gain');
subplot(2,2,3);
scatter(DoGtest.rfs,DoGtest.normamp)
xlabel('rfs');ylabel('namp');
subplot(2,2,4);
scatter(DoGtest.rfs,DoGtest.sdratio)
xlabel('rfs');ylabel('sdratio');
close(ftemp)

% calculate the actual pRF shape
x=linspace(-20,20,1000);
y = (...
    (1./(2*pi*DoGtest.rfs).*exp(-(x).^2./(2*DoGtest.rfs.^2))) - ...
    DoGtest.normamp.*(1./(2*pi*(DoGtest.rfs.*DoGtest.sdratio)).*...
    exp(-(x).^2./(2*(DoGtest.rfs.*DoGtest.sdratio).^2))) ...
    );
f_dog = figure;
plot(x,y);

% mean
ym = nanmean(y,1); 
hold on; plot(x,ym,'k','Linewidth',5)
close(f_dog)


fprintf('Handling DoG pRF sizes MUA \n')
R2th=70;
DoGtest = tMUA(...
    strcmp(tMUA.Model,'dog_ephys_cv1') & ...
    tMUA.Area == 1 & ...
    tMUA.R2 > R2th ,:);

ftemp = figure;
subplot(2,2,1);
scatter(DoGtest.ecc,DoGtest.rfs);
xlabel('ecc');ylabel('rfs');
subplot(2,2,2);
scatter(DoGtest.rfs,DoGtest.gain)
xlabel('rfs');ylabel('gain');
subplot(2,2,3);
scatter(DoGtest.rfs,DoGtest.normamp)
xlabel('rfs');ylabel('namp');
subplot(2,2,4);
scatter(DoGtest.rfs,DoGtest.sdratio)
xlabel('rfs');ylabel('sdratio');
close(ftemp)

% calculate the actual pRF shape
x=linspace(-20,20,1000);
y = (...
    (1./(2*pi*DoGtest.rfs).*exp(-(x).^2./(2*DoGtest.rfs.^2))) - ...
    DoGtest.normamp.*(1./(2*pi*(DoGtest.rfs.*DoGtest.sdratio)).*...
    exp(-(x).^2./(2*(DoGtest.rfs.*DoGtest.sdratio).^2))) ...
    );
f_dog = figure;
plot(x,y./max(y,2));

% mean
ym = nanmean(y,1); 
hold on; plot(x,ym,'k','Linewidth',5)
close(f_dog)





%% ========================================================================
%  FMRI PRF ANALYSIS ------------------------------------------------------
%  ========================================================================

%% Number of significant voxels per ROI (split by monkey and model) =======
RTHRES = 5;

ff = figure;
set(ff,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
set(ff,'Position',[10 10 1600 1000]);

paneln=1; 
for s = 1:length(SUBS)
    fprintf(['SUB: ' SUBS{s} '\n'])
    for rowmod=1:4
            subplot(4,2,paneln); hold on;
            nvox=[];
            for r=1:length(roi)
                nvox = [nvox; sum(...
                    T(modidx.(MRI_MODELS{rowmod,1})).mod.R2 > RTHRES & ...
                    strcmp(T(modidx.(MRI_MODELS{rowmod,1})).mod.Monkey,SUBS{s}) & ...
                    ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois) ) ) ...
                     sum(strcmp(T(modidx.(MRI_MODELS{rowmod,1})).mod.Monkey,SUBS{s}) & ...
                    ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois) ) ) ];
            end
            for xval=1:size(nvox,1)
                bar(xval,nvox(xval,1));
            end
            title(['SUB ' SUBS{s} ', ' MMS{rowmod,1} ' nVox with R2 > ' ...
                num2str(RTHRES)],'interpreter','none');
            xlabel('ROI'); ylabel('nVox','interpreter','none');
            set(gca, 'Box','off','yscale','log');
            set(gca,'xtick',1:length(nvox),'xticklabels',roilabels,'TickDir','out');
            xtickangle(45)
            paneln=paneln+1;
    end
end
saveas(ff,fullfile(pngfld, 'MRI_nVoxSign_ROI.png'));
saveas(ff,fullfile(svgfld, 'MRI_nVoxSign_ROI.svg'));
if CloseFigs; close(ff); end
    
ff2 = figure;
set(ff2,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
set(ff2,'Position',[10 10 1600 1000]);

paneln=1; 
for s = 1: length(SUBS)
    fprintf(['SUB: ' SUBS{s} '\n'])
    for rowmod=1:4
            subplot(4,2,paneln); hold on;
            nvox=[];
            for r=1:length(roi)
                nvox = [nvox; sum(...
                    T(modidx.(MRI_MODELS{rowmod,1})).mod.R2 > RTHRES & ...
                    strcmp(T(modidx.(MRI_MODELS{rowmod,1})).mod.Monkey,SUBS{s}) & ...
                    ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois) ) ) ...
                     sum(strcmp(T(modidx.(MRI_MODELS{rowmod,1})).mod.Monkey,SUBS{s}) & ...
                    ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois) ) ) ];
            end
            for xval=1:size(nvox,1)
                bar(xval,nvox(xval,1)./nvox(xval,2));
            end
            title(['SUB ' SUBS{s} ', ' MMS{rowmod,1} ' proportion Vox with R2 > ' ...
                num2str(RTHRES)],'interpreter','none');
            xlabel('ROI'); ylabel('Proportion Vox','interpreter','none');
            set(gca, 'Box','off','ylim',[0 .6]);
            set(gca,'xtick',1:length(nvox),'xticklabels',roilabels,'TickDir','out');
            xtickangle(45)
            paneln=paneln+1;
    end
end

if SaveFigs; saveas(ff2,fullfile(pngfld, 'MRI_propVoxSign_ROI.png')); end
if SaveFigs; saveas(ff2,fullfile(svgfld, 'MRI_propVoxSign_ROI.svg')); end
if CloseFigs; close(ff2); end

%% MRI scatter plots & differences R2 per ROI =============================
RTHRES = 0; % only include when one of the models fits above threshold

for sidx = 0:2 % both monkeys and individuals   
    % scatter plots ---
    % Different ROIS in different colors
    f=figure;
    set(f,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
    set(f,'Position',[10 10 1300 600]);
    
    % Tell us what's happening
    if sidx
        fprintf(['Monkey ' SUBS{sidx} '\n']);
        marker = [' ' SUBS{sidx}];
    else
        fprintf(['Both monkeys\n']);
        marker = ' both';
    end
    
    % plot scatter per area
    figure(f); paneln=1;
    for rowmod=1:4
        for colmod=rowmod+1:4
            subplot(2,3,paneln); hold on;
            PerROI(sidx+1,paneln).Models = {MRI_MODELS{rowmod,1},MRI_MODELS{colmod,1}};
            PerROI(sidx+1,paneln).mR2 = []; PerROI(sidx+1,paneln).seR2 = [];
            for r=1:length(roi)
                s_R2 = T(modidx.(MRI_MODELS{rowmod,1})).mod.R2 >= RTHRES  | ...
                    T(modidx.(MRI_MODELS{colmod,1})).mod.R2 >= RTHRES;
                if sidx
                    s_Monkey = strcmp(T(modidx.(MRI_MODELS{rowmod,1})).mod.Monkey,SUBS{sidx});
                    SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois)) & s_Monkey;
                else
                    SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois));
                end
                scatter(...
                    T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS),...
                    T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS),...
                    100,'Marker','.');
                PerROI(sidx+1,paneln).roi(r).sub = marker(2:end);
                PerROI(sidx+1,paneln).roi(r).name = roilabels(r);
                PerROI(sidx+1,paneln).roi(r).R2 = [...
                    T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS),...
                    T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS)];
                PerROI(sidx+1,paneln).mR2 = [PerROI(sidx+1,paneln).mR2; ...
                    nanmean(T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS)),...
                    nanmean(T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS))];
                PerROI(sidx+1,paneln).seR2 = [PerROI(sidx+1,paneln).seR2; ...
                    nanstd(T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS))./...
                    sqrt(length(T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS))),...
                    nanstd(T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS))./...
                    sqrt(length(T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS)))];
                % stats per ROI
                [PerROI(sidx+1,paneln).roi(r).p,...
                    PerROI(sidx+1,paneln).roi(r).h,...
                    PerROI(sidx+1,paneln).roi(r).stats] = signrank(...
                    PerROI(sidx+1,paneln).roi(r).R2(:,1),...
                    PerROI(sidx+1,paneln).roi(r).R2(:,2));
                PerROI(sidx+1,paneln).roi(r).dir = ...
                    diff(PerROI(sidx+1,paneln).mR2(r,:))./...
                    abs(diff(PerROI(sidx+1,paneln).mR2(r,:)));
            end
            plot([-2 100],[-2 100],'k','Linewidth',1);
            title([MMS{rowmod,1} ' vs ' MMS{colmod,1}],'interpreter','none');
            xlabel(MMS{rowmod,1},'interpreter','none');
            ylabel(MMS{colmod,1},'interpreter','none');
            set(gca, 'Box','off', 'xlim', [-2 100], 'ylim',[-2 100]);
            paneln=paneln+1;
        end
    end
    suptitle(['SUBJECT' marker])
    if SaveFigs; saveas(f,fullfile(pngfld, ['MRI_ModelComparison_ROI_R2' marker '.png'])); end
    if SaveFigs; saveas(f,fullfile(svgfld, ['MRI_ModelComparison_ROI_R2' marker '.svg'])); end
    if CloseFigs; close(f); end
end

% Report stats per ROI ----
fprintf('Compare CSS against P-LIN per ROI -----\n');
nBetter = sum([PerROI(1,2).roi(:).p] < 0.05 & [PerROI(1,2).roi(:).dir] > 0);
nWorse = sum([PerROI(1,2).roi(:).p] < 0.05 & [PerROI(1,2).roi(:).dir] < 0);
nSimilar = sum([PerROI(1,2).roi(:).p] > 0.05);
nROI = length(PerROI(1,2).roi);
fprintf('- Both together -\n');
fprintf([num2str(nBetter) '/' num2str(nROI) ' CSS better\n']);
fprintf([num2str(nWorse) '/' num2str(nROI) ' CSS worse\n']);
fprintf([num2str(nSimilar) '/' num2str(nROI) ' no difference\n']);

nBetter = sum([PerROI(2,2).roi(:).p] < 0.05 & [PerROI(2,2).roi(:).dir] > 0);
nWorse = sum([PerROI(2,2).roi(:).p] < 0.05 & [PerROI(2,2).roi(:).dir] < 0);
nSimilar = sum([PerROI(2,2).roi(:).p] > 0.05);
nROI = length(PerROI(2,2).roi);
fprintf(['- ' SUBS{1} ' -\n']);
fprintf([num2str(nBetter) '/' num2str(nROI) ' CSS better\n']);
fprintf([num2str(nWorse) '/' num2str(nROI) ' CSS worse\n']);
fprintf([num2str(nSimilar) '/' num2str(nROI) ' no difference\n']);

nBetter = sum([PerROI(3,2).roi(:).p] < 0.05 & [PerROI(3,2).roi(:).dir] > 0);
nWorse = sum([PerROI(3,2).roi(:).p] < 0.05 & [PerROI(3,2).roi(:).dir] < 0);
nSimilar = sum([PerROI(3,2).roi(:).p] > 0.05);
nROI = length(PerROI(3,2).roi);
fprintf(['- ' SUBS{2} ' -\n']);
fprintf([num2str(nBetter) '/' num2str(nROI) ' CSS better\n']);
fprintf([num2str(nWorse) '/' num2str(nROI) ' CSS worse\n']);
fprintf([num2str(nSimilar) '/' num2str(nROI) ' no difference\n']);

%% MRI scatter plots & differences R2 all voxels ==========================
RTHRES = 0; % only include when one of the models fits above threshold

for sidx = 0:2 % both monkeys and individuals   

    % All ROIS together in black
    fsc=figure;
    set(fsc,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
    set(fsc,'Position',[10 10 1300 600]);
    
    % Tell us what's happening
    if sidx
        fprintf(['Monkey ' SUBS{sidx} '\n']);
        marker = [' ' SUBS{sidx}];
    else
        fprintf(['Both monkeys\n']);
        marker = ' both';
    end
    
    % plot scatter over all voxels
    figure(fsc); paneln=1;
    for rowmod=1:4
        for colmod=rowmod+1:4
            subplot(2,3,paneln); hold on;
            XY=[];
            for r=1:length(roi) 
                s_R2 = T(modidx.(MRI_MODELS{rowmod,1})).mod.R2 >= RTHRES  | ...
                    T(modidx.(MRI_MODELS{colmod,1})).mod.R2 >= RTHRES;
                if sidx
                    s_Monkey = strcmp(T(modidx.(MRI_MODELS{rowmod,1})).mod.Monkey,SUBS{sidx});
                    SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois)) & s_Monkey;
                else
                    SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois));
                end
                XY=[XY; [T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS) ...
                    T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS)] ];
            end
            binscatter(XY(:,1),XY(:,2),100); colorbar; 
            set(gca,'ColorScale','log');
            colormap(inferno)
            plot([-2 100],[-2 100],'k','Linewidth',1);
            title([MMS{rowmod,1} ' vs ' MMS{colmod,1}],'interpreter','none');
            xlabel(MMS{rowmod,1},'interpreter','none');
            ylabel(MMS{colmod,1},'interpreter','none');
            set(gca, 'Box','off', 'xlim', [-2 100], 'ylim',[-2 100]);
            caxis([1 1e4])  
            paneln=paneln+1;
        end
    end
    suptitle(['SUBJECT' marker])
    if SaveFigs; saveas(fsc,fullfile(pngfld, ['MRI_ModelComparison_R2' marker '.png'])); end
    if SaveFigs; saveas(fsc,fullfile(svgfld, ['MRI_ModelComparison_R2' marker '.svg'])); end
    if CloseFigs; close(fsc); end

end

% stats
MR2 = [T(2).mod.R2 T(4).mod.R2 T(7).mod.R2 T(8).mod.R2];
sel=logical(sum(MR2>0,2));
[p,tbl,stats] = kruskalwallis(MR2(sel,:),{'css','dog','p-lin','u-lin'});
[c,m,h,gnames] = multcompare(stats);
for i=1:size(c,1)
    fprintf([gnames{c(i,1)} ' vs ' gnames{c(i,2)} ...
        ', p = ' num2str(c(i,6))  '\n'])
end



% diff distributions plots -----
for sidx = 0:2 % both monkeys and individuals
    % Tell us what's happening
    if sidx
        fprintf(['Monkey ' SUBS{sidx} '\n']);
        marker = [' ' SUBS{sidx}];
    else
        fprintf(['Both monkeys\n']);
        marker = ' both';
    end
    
    
    idx=1;
    for rowmod=1:4
        for colmod=rowmod+1:4
            diffmat{sidx+1,idx}=[];
            for r=1:length(roi)               
                if sidx
                    s_Monkey = strcmp(T(modidx.(MRI_MODELS{rowmod,1})).mod.Monkey,SUBS{sidx});
                    SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois)) & s_Monkey;
                else
                    SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                        ck_GetROIidx(roilabels(r),rois));
                end
                [n,x] = hist(...
                    T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS)-...
                    T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS), 100);
                f = n./sum(SSS);
                
                m = mean(T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS)-...
                    T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS));
                sd = std(T(modidx.(MRI_MODELS{colmod,1})).mod.R2(SSS)-...
                    T(modidx.(MRI_MODELS{rowmod,1})).mod.R2(SSS));
                se = sd ./ sqrt(sum(SSS));
                diffmat{sidx+1,idx} = [diffmat{sidx+1,idx}; m sd se];
            end
            idx=idx+1;
        end
    end
end

    
f2=figure;
set(f2,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
set(f2,'Position',[100 100 1800 1000]);

cc=1;
for rowmod=1:4
    for colmod=rowmod+1:4
        subplot(3,2,cc); hold on;
        for xval=1:length(diffmat{1,cc})
            bar(xval,diffmat{1,cc}(xval,1));
        end
        for xval=1:length(diffmat{1,cc})
            errorbar(xval,diffmat{1,cc}(xval,1),diffmat{cc}(xval,3),...
                'k-','Linestyle','none');
        end
        set(gca,'xtick',1:length(diffmat{cc}),...
            'xticklabels',roilabels,'ylim',[-2 2.5],'TickDir','out');
        xlabel('ROI'); ylabel('Diff R2');
        title([MMS{colmod,1} ' - ' MMS{rowmod,1}],...
            'interpreter','none');
        xtickangle(45)
        cc=cc+1;
    end
end
if SaveFigs; saveas(f2,fullfile(figfld, 'MRI_ModelComparison_ROI_R2.png')); end
if CloseFigs; close(f2); end

%% only CSS and DoG vs P-LIN

f3=figure;
set(f3,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
set(f3,'Position',[100 100 1200 1000]);

%subplot(3,1,1); hold on;
subplot(2,2,1); hold on;
for xval=1:length(diffmat{2,2})
    bar(xval,diffmat{2,2}(xval,1));
end
for xval=1:length(diffmat{2,2})
    errorbar(xval,diffmat{2,2}(xval,1),diffmat{2,2}(xval,3),...
        'k-','Linestyle','none');
end
set(gca,'xtick',1:length(diffmat{2,2}),...
    'xticklabels',roilabels,'ylim',[-0.1 2.5],'TickDir','out');
xlabel('ROI'); ylabel('Diff R2');
title('M1: CSS - P-LIN','interpreter','none');
xtickangle(45)

%subplot(3,1,2); hold on;
subplot(2,2,2); hold on;

for xval=1:length(diffmat{2,3})
    bar(xval,diffmat{2,3}(xval,1));
end
for xval=1:length(diffmat{2,3})
    errorbar(xval,diffmat{2,3}(xval,1),diffmat{2,3}(xval,3),...
        'k-','Linestyle','none');
end
set(gca,'xtick',1:length(diffmat{2,3}),...
    'xticklabels',roilabels,'ylim',[-0.1 2.5],'TickDir','out');
xlabel('ROI'); ylabel('Diff R2');
title('M1: DoG - P-LIN','interpreter','none');
xtickangle(45)


subplot(2,2,3); hold on;
for xval=1:length(diffmat{3,2})
    bar(xval,diffmat{3,2}(xval,1));
end
for xval=1:length(diffmat{3,2})
    errorbar(xval,diffmat{3,2}(xval,1),diffmat{3,2}(xval,3),...
        'k-','Linestyle','none');
end
set(gca,'xtick',1:length(diffmat{3,2}),...
    'xticklabels',roilabels,'ylim',[-0.1 2.5],'TickDir','out');
xlabel('ROI'); ylabel('Diff R2');
title('M2: CSS - P-LIN','interpreter','none');
xtickangle(45)

%subplot(3,1,2); hold on;
subplot(2,2,4); hold on;
for xval=1:length(diffmat{3,3})
    bar(xval,diffmat{3,3}(xval,1));
end
for xval=1:length(diffmat{3,3})
    errorbar(xval,diffmat{3,3}(xval,1),diffmat{3,3}(xval,3),...
        'k-','Linestyle','none');
end
set(gca,'xtick',1:length(diffmat{3,3}),...
    'xticklabels',roilabels,'ylim',[-0.1 2.5],'TickDir','out');
xlabel('ROI'); ylabel('Diff R2');
title('M2: DoG - P-LIN','interpreter','none');
xtickangle(45)





% subplot(3,1,3); hold on;
% for xval=1:length(diffmat{6})
%     bar(xval,diffmat{6}(xval,1));
% end
% for xval=1:length(diffmat{6})
%     errorbar(xval,diffmat{6}(xval,1),diffmat{6}(xval,3),...
%         'k-','Linestyle','none');
% end
% set(gca,'xtick',1:length(diffmat{6}),...
%     'xticklabels',roilabels,'ylim',[-2.1 2.1],'TickDir','out');
% xlabel('ROI'); ylabel('Diff R2');
% title('DoG - CSS','interpreter','none');
% xtickangle(45)

if SaveFigs; saveas(f3,fullfile(figfld, 'MRI_SelModelComparison_ROI_R2.png')); end
if CloseFigs; close(f3); end

%% Good DoG and NegGain Fits: Characterize ================================
R2th = 5; % minimum R2
R2enh = 5; % R2 improvement

% DoG = tMRI(...
%     strcmp(tMRI.Model,'doghrf_cv1_mhrf'),:);
% lin_n = tMRI(...
%     strcmp(tMRI.Model,'linhrf_cv1_mhrf_neggain'),:);
% lin = tMRI(...
%     strcmp(tMRI.Model,'linhrf_cv1_mhrf'),:);
% 
% f_neg1 = figure;
% set(f_neg1,'Position',[10 10 1200 1000]);
% vox_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh;
% subplot(2,2,1);scatter(DoG.X(vox_sel),DoG.Y(vox_sel),...
%     'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
% set(gca,'xaxislocation','origin','yaxislocation','origin',...
%     'xlim',[-8 8],'ylim',[-8 8]);
% title('Locations of pRF with good DoG fits & bad LIN fits')
% xlabel('X deg');ylabel('Y deg');
% 
% subplot(2,2,2);histogram(DoG.ecc(vox_sel),0:0.1:5,...
%     'FaceColor','k','FaceAlpha',0.5);
% title('ECC of pRF with good DoG fits & bad LIN fits')
% ylabel('nvoxels'); xlabel('Ecc');
% 
% vox_sel = lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh;
% subplot(2,2,3);scatter(lin_n.X(vox_sel),lin_n.Y(vox_sel),...
%     'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
% set(gca,'xaxislocation','origin','yaxislocation','origin',...
%     'xlim',[-8 8],'ylim',[-8 8]);
% title('Locations of pRF with good LIN-N fits & bad LIN fits')
% xlabel('X deg');ylabel('Y deg');
% 
% subplot(2,2,4);histogram(lin_n.ecc(vox_sel),0:0.1:5,...
%     'FaceColor','k','FaceAlpha',0.5);
% title('ECC of pRF with good LIN-N fits & bad LIN fits')
% ylabel('nvoxels'); xlabel('Ecc');
% 
% if SaveFigs; saveas(f_neg1,fullfile(figfld, 'MRI_NEG-PRF1.png')); end
% if CloseFigs; close(f_neg1); end

% ---------

% f_neg2 = figure;
% set(f_neg2,'Position',[10 10 1300 1600]);
% vox_sel = lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh & ...
%     (lin_n.ecc<15 & lin.ecc<15 & DoG.ecc<15) & ....
%     (lin_n.rfs<10 & lin.rfs<10 & DoG.rfs<10);
% 
% subplot(3,2,1); hold on;
% plot([0 15],[0 15],'r');
% scatter(lin_n.ecc(vox_sel),lin.ecc(vox_sel),...
%     'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
% xlabel('Ecc. LINEAR POSNEG')
% ylabel('Ecc. LINEAR POS')
% set(gca,'xlim',[0 15],'ylim',[0 15]);
% yy=get(gca,'ylim');
% title('Eccentricity')
% 
% subplot(3,2,2); hold on;
% histogram(lin_n.gain(vox_sel),-10:0.1:10,'FaceColor','k','FaceAlpha',0.5);
% xlabel('gain LIN-POSNEG');ylabel('nvoxels');
% set(gca,'xlim',[-3 1]);
% MM=median(lin_n.gain(vox_sel));
% yy=get(gca,'ylim');
% plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
% set(gca,'ylim',[0 yy(2)+30]);
% title('Gain')
% 
% subplot(3,2,3); hold on;
% bb = [lin_n.ecc(vox_sel) lin.ecc(vox_sel)];
% plot([1 2],bb,'LineWidth',0.1,'Color',[.8 .8 .8])
% plot([1 2],mean(bb),'k','Linewidth',5)
% set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
%     'ylim',[0 20],'xlim',[0.8 2.2])
% ylabel('Eccentricity');
% title('Ecc Diff')
% 
% subplot(3,2,4); hold on;
% histogram(lin.ecc(vox_sel)-lin_n.ecc(vox_sel),-10:0.5:10,...
%     'FaceColor','k','FaceAlpha',0.5);
% xlabel('Ecc. Diff (POS-POSNEG)');ylabel('nvoxels');
% MM=median(bb(:,2)-bb(:,1));
% yy=get(gca,'ylim');
% plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
% set(gca,'ylim',[0 yy(2)+30]);
% title('Ecc Diff')
% 
% subplot(3,2,5); hold on;
% bb = [lin_n.rfs(vox_sel) lin.rfs(vox_sel)];
% plot([1 2],bb)
% plot([1 2],mean(bb),'k','Linewidth',5)
% set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
%     'ylim',[0 6],'xlim',[0.8 2.2])
% ylabel('Size');
% title('Size Diff')
% 
% subplot(3,2,6); hold on;
% histogram(lin.rfs(vox_sel)-lin_n.rfs(vox_sel),-10:0.5:10,...
%     'FaceColor','k','FaceAlpha',0.5);
% xlabel('Size Diff (POS-POSNEG)');ylabel('nvoxels');
% MM=median(bb(:,2)-bb(:,1));
% yy=get(gca,'ylim');
% plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
% set(gca,'ylim',[0 yy(2)+30]);
% set(gca,'xlim',[-5 5]);
% title('Size Diff')
% 
% sgtitle('pRFs POS LINEAR vs POSNEG LINEAR model')
% 
% if SaveFigs; saveas(f_neg2,fullfile(figfld, 'MRI_NEG-PRF2.png')); end
% if CloseFigs; close(f_neg2); end

% ---------

% f_neg3 = figure;
% set(f_neg3,'Position',[10 10 1300 1100]);
% vox_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh;
% 
% subplot(2,2,1); hold on;
% plot([0 15],[0 15],'r');
% scatter(DoG.ecc(vox_sel),lin.ecc(vox_sel),...
%     'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
% xlabel('Ecc. DOG')
% ylabel('Ecc. LINEAR POS')
% set(gca,'xlim',[0 15],'ylim',[0 15]);
% yy=get(gca,'ylim');
% plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
% title('Eccentricity')
% 
% subplot(2,2,2); hold on;
% histogram(DoG.normamp(vox_sel),-10:0.1:10,'FaceColor','k','FaceAlpha',0.5);
% xlabel('INH nAMP');ylabel('nvoxels');
% set(gca,'xlim',[0 3]);
% MM=median(DoG.normamp(vox_sel));
% yy=get(gca,'ylim');
% plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
% set(gca,'ylim',[0 yy(2)+30]);
% title('NORMAMP')
% 
% subplot(2,2,3); hold on;
% bb = [DoG.ecc(vox_sel) lin.ecc(vox_sel)];
% plot([1 2],bb)
% plot([1 2],mean(bb),'k','Linewidth',5)
% set(gca,'xtick',1:2,'xticklabels',{'DoG','LIN'},...
%     'ylim',[0 20],'xlim',[0.8 2.2])
% ylabel('Eccentricity');
% title('Ecc Diff')
% 
% subplot(2,2,4); hold on;
% histogram(lin.ecc(vox_sel)-DoG.ecc(vox_sel),-10:0.5:10,...
%     'FaceColor','k','FaceAlpha',0.5);
% xlabel('Ecc. Diff (POS-DoG)');ylabel('nvoxels');
% MM=median(bb(:,2)-bb(:,1));
% yy=get(gca,'ylim');
% plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
% set(gca,'ylim',[0 yy(2)+30]);
% title('Ecc Diff')
% 
% sgtitle('pRFs POS LINEAR vs DoG model')
% 
% if SaveFigs; saveas(f_neg3,fullfile(figfld, 'MRI_NEG-PRF3.png')); end
% if CloseFigs; close(f_neg3); end

% ---------

DoG = tMRI(...
    strcmp(tMRI.Model,'doghrf_cv1_mhrf'),:);
lin_n = tMRI(...
    strcmp(tMRI.Model,'linhrf_cv1_mhrf_neggain'),:);
lin = tMRI(...
    strcmp(tMRI.Model,'linhrf_cv1_mhrf'),:);

f_neg2 = figure;
set(f_neg2,'Position',[10 10 1200 1000]);

vox_sel = lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh & ...
    (lin_n.ecc<25 & lin.ecc<25 & DoG.ecc<25) & ....
    (lin_n.rfs<20 & lin.rfs<20 & DoG.rfs<20);

fprintf('-------------\n');

subplot(2,3,1); hold on;
histogram(lin_n.gain(vox_sel),-100:0.1:100,...
    'Normalization','probability','FaceColor','k','FaceAlpha',0.75);
xlabel('gain LIN-POSNEG');ylabel('nvoxels');
set(gca,'xlim',[-2.6 2.6]);
MM=median(lin_n.gain(vox_sel));
yy=get(gca,'ylim');
plot([MM MM], [0 1],'k','Linewidth',2)
set(gca,'ylim',[0 0.3],'TickDir','out');
title('Gain SELECTED voxels')
fprintf(['MEDIAN GAIN SEL: ' num2str(MM) '\n'])

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(lin_n.gain(vox_sel),0,'tail','left');
fprintf(['Gain SEL < 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

% all vox
subplot(2,3,4); hold on;
histogram(lin_n.gain,-100:0.1:100,...
    'Normalization','probability','FaceColor','k','FaceAlpha',0.75);
xlabel('gain LIN-POSNEG');ylabel('nvoxels');
set(gca,'xlim',[-2.6 2.6]);
MM=median(lin_n.gain);
yy=get(gca,'ylim');
plot([MM MM], [0 1],'k','Linewidth',2)
set(gca,'ylim',[0 0.3],'TickDir','out');
title('Gain ALL voxels')
fprintf(['MEDIAN GAIN ALL: ' num2str(MM) '\n'])

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(lin_n.gain,0,'tail','right');
fprintf(['Gain ALL < 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf('-------------\n');

subplot(2,3,2); hold on;
bb = [lin_n.ecc(vox_sel) lin.ecc(vox_sel)];
plot([1 2],mean(bb),'o','Linewidth',2)
errorbar([1 2],mean(bb),std(bb),'k','Linewidth',2)
plot([1 2],mean(bb),'o','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k')
set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
    'ylim',[0 15],'xlim',[0.8 2.2],'TickDir','out')
ylabel('Eccentricity');
title('Ecc Diff')

subplot(2,3,3); hold on;
histogram(lin.ecc(vox_sel)-lin_n.ecc(vox_sel),-10:0.5:10,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Ecc. Diff (POS-POSNEG)');ylabel('nvoxels');
MM=median(bb(:,2)-bb(:,1));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+100],'TickDir','out');
title('Ecc Diff')
fprintf(['MEDIAN ECC DIFF: ' num2str(MM) '\n'])

bbecc=bb;

subplot(2,3,5); hold on;
bb = [lin_n.rfs(vox_sel) lin.rfs(vox_sel)];
plot([1 2],mean(bb),'o','Linewidth',2)
errorbar([1 2],mean(bb),std(bb),'k','Linewidth',2)
plot([1 2],mean(bb),'o','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k')
set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
    'ylim',[0 5.1],'xlim',[0.8 2.2],'TickDir','out')
ylabel('Size');
title('Size Diff')

subplot(2,3,6); hold on;
histogram(lin.rfs(vox_sel)-lin_n.rfs(vox_sel),-10:0.25:10,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Size Diff (POS-POSNEG)');ylabel('nvoxels');
MM=median(bb(:,2)-bb(:,1));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+110]);
set(gca,'xlim',[-5.1 5.1],'TickDir','out');
title('Size Diff')
fprintf(['MEDIAN SZ DIFF: ' num2str(MM) '\n'])

bbsz=bb;
sgtitle('pRFs POS LINEAR vs POSNEG LINEAR model')

if SaveFigs; saveas(f_neg2,fullfile(pngfld, 'MRI_NEG-PRF2.png')); end
if SaveFigs; saveas(f_neg2,fullfile(svgfld, 'MRI_NEG-PRF2.svg')); end
if CloseFigs; close(f_neg2); end

% Wilcoxon 1-tailed < 0
[p,h,stats] = signrank(bbecc(:,1),bbecc(:,2));
fprintf(['Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

%% ---------

f_neg3 = figure;
set(f_neg3,'Position',[10 10 1300 1100]);
vox_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh;

subplot(1,3,1); hold on;
histogram(DoG.normamp(vox_sel),-10:0.1:10,'FaceColor','k','FaceAlpha',0.5);
xlabel('INH nAMP');ylabel('nvoxels');
set(gca,'xlim',[-0.2 3.1]);
MM=median(DoG.normamp(vox_sel));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+30],'TickDir','out');
title('NORMAMP')
fprintf(['MEDIAN NAMP: ' num2str(MM) '\n'])


% Wilcoxon 1-tailed > 1
[p,h,stats] = signrank(DoG.normamp(vox_sel),1,'tail','right');
fprintf(['nAmp > 1: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['Median nAMP: ' num2str(median(DoG.normamp(vox_sel))) ', IQR: '...
    num2str(iqr(DoG.normamp(vox_sel))) '\n'])

subplot(1,3,2); hold on;
bb = [DoG.ecc(vox_sel) lin.ecc(vox_sel)];
bbecc=bb;
plot([1 2],mean(bb),'o','Linewidth',2)
errorbar([1 2],mean(bb),std(bb),'k','Linewidth',2)
plot([1 2],mean(bb),'o','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k')
set(gca,'xtick',1:2,'xticklabels',{'DoG','LIN'},...
    'ylim',[0 20],'xlim',[0.8 2.2],'TickDir','out')
ylabel('Eccentricity');
title('Ecc Diff')

subplot(1,3,3); hold on;
histogram(lin.ecc(vox_sel)-DoG.ecc(vox_sel),-10:0.5:10,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Ecc. Diff (POS-DoG)');ylabel('nvoxels');
MM=median(bb(:,2)-bb(:,1));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+110],'TickDir','out');
title('Ecc Diff')
fprintf(['MEDIAN ECC DIFF: ' num2str(MM) '\n'])

sgtitle('pRFs POS LINEAR vs DoG model')

%if SaveFigs; saveas(f_neg3,fullfile(figfld, 'MRI_NEG-PRF3.png')); end
%if CloseFigs; close(f_neg3); end

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(bbecc(:,1),bbecc(:,2));
fprintf(['Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

%% Value of exponential parameter for CSS across ROIs =====================
RTHRES = 5;
fexp = figure;
set(fexp,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
set(fexp,'Position',[10 10 1600 300]);
hold on;

for s = 1:length(SUBS)
    fprintf(['SUB: ' SUBS{s} '\n'])
    exptv(s).roi={};
    for r=1:length(roi)
        exptvals = T(modidx.csshrf_cv1_mhrf).mod.expt(...
            T(modidx.csshrf_cv1_mhrf).mod.R2 > RTHRES & ...
            strcmp(T(modidx.csshrf_cv1_mhrf).mod.Monkey,SUBS{s}) & ...
            ismember(T(modidx.csshrf_cv1_mhrf).mod.ROI,...
            ck_GetROIidx(roilabels(r),rois) ) );
        exptv(s).roi{r}=exptvals;
    end
end

all_expt=[];
for xval=1:length(roi)
    if size([exptv(1).roi{xval};exptv(2).roi{xval}],1) > 4
        bar(xval,mean([exptv(1).roi{xval};exptv(2).roi{xval}]));
        all_expt = [all_expt; exptv(1).roi{xval}; exptv(2).roi{xval}];
    else
        bar(xval,0);
        all_expt = [all_expt; exptv(1).roi{xval}; exptv(2).roi{xval}];
    end
end

for xval=1:length(roi)
    if size([exptv(1).roi{xval};exptv(2).roi{xval}],1) > 4
        errorbar(xval,...
            mean([exptv(1).roi{xval};exptv(2).roi{xval}]),...
            std([exptv(1).roi{xval};exptv(2).roi{xval}]),...
            'k','Linestyle','none');
        %     errorbar(xval,...
        %         mean([exptv(1).roi{xval};exptv(2).roi{xval}]),...
        %         std([exptv(1).roi{xval};exptv(2).roi{xval}])./...
        %         sqrt(length([exptv(1).roi{xval};exptv(2).roi{xval}])),...
        %         'k','Linestyle','none');
    end
end
title(['Avg EXPT parameter, R2 > ' ...
    num2str(RTHRES)],'interpreter','none');
%xlabel('ROI'); 
ylabel('EXPT PM','interpreter','none');
set(gca,'xtick',1:length(roi),'xticklabels',roilabels,'TickDir','out');
xtickangle(45)
if SaveFigs; saveas(fexp,fullfile(figfld, 'MRI_CSS_EXPT-PM_ROI.png')); end
if CloseFigs; close(fexp); end

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(all_expt,1,'tail','left');
fprintf(['ALL VOXELS - Wilcoxon EXPT < 1 : z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

% [p,h,stats] = signrank([exptv(1).roi{1};exptv(2).roi{1}],1,'tail','left');
% fprintf(['\nV1 VOXELS - Wilcoxon EXPT < 1 : z = ' ...
%     num2str(stats.zval) ', p = ' num2str(p) '\n\n']);


fprintf('\n------\nSTATS\n------\n')
for xval=1:length(roi)
    if size([exptv(1).roi{xval};exptv(2).roi{xval}],1) > 4
        [p,h,stats] = signrank([exptv(1).roi{xval};exptv(2).roi{xval}],1,'tail','left');
        if isfield(stats,'zval')  
            fprintf([roilabels{xval} ' VOXELS - Wilcoxon EXPT < 1 : z = ' ...
                num2str(stats.zval) ', p = ' num2str(p) '\n']);
        end
    else
        fprintf([roilabels{xval} ': Not enough values for statistics\n'])
    end
end

%% MRI scatter plot HRF & differences =====================================
RTHRES = 0;

% f3=figure;
% set(f3,'Position',[100 100 2000 1200]);
% set(f3,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
% s_R2 = T(modidx.linhrf_cv1_mhrf).mod.R2>RTHRES;
% 
% for mm = 1:size(MRI_MODELS,1)
%     subplot(4,3,(mm-1)*3 +1); hold on;
%     plot([0 100],[0 100],'k','LineWidth',2);
%     for r=1:length(roi)
%         SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
%                 ck_GetROIidx(roilabels(r),rois) );
%         scatter(T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS),...
%             T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS),100,'Marker','.');
%     end
%     title(['mHRF vs cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
%     xlabel('Monkey HRF R2'); ylabel('Canonical HRF R2');
%     set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
% 
%     diffmat2{1}=[];
%     for r=1:length(roi)
%         SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
%                 ck_GetROIidx(roilabels(r),rois) );
%         [n,x] = hist(...
%             T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS)-...
%             T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS),100);
%         f = n./sum(SSS);
%     
%         m = mean(T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS)-...
%             T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS));
%         sd = std(T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS)-...
%             T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS));
%         nvox = sum(SSS);
%         se = sd ./ sqrt(nvox);
%         diffmat2{1} = [diffmat2{1}; m sd se nvox];
%     end
% 
%     subplot(4,3,(mm-1)*3 +2); hold on
%     for xval=1:length(diffmat2{1})
%         bar(xval,diffmat2{1}(xval,1));
%     end
%     for xval=1:length(diffmat2{1})
%         errorbar(xval,diffmat2{1}(xval,1),diffmat2{1}(xval,3),...
%         'k-','Linestyle','none')
%     end
%     % mean (taking nvox into account)
%     mAll = sum((diffmat2{1}(:,1).*diffmat2{1}(:,4)))./...
%         sum(diffmat2{1}(:,4));    
%     text(0.5, -0.4,num2str(mAll))
%     set(gca,'xticklabels',[],'ylim',[-0.65 1]);
%     xlabel('ROI'); ylabel('Diff R2');
%     title(['mHRF - cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
%     %legend(roilabels,'interpreter','none','Location','NorthEast');
%     set(gca,'xtick',1:length(diffmat2{1}),'xticklabels',roilabels);
%     xtickangle(45)
%     
%     subplot(4,3,(mm-1)*3 +3); hold on
%     for xval=1:length(diffmat2{1})
%         bar(xval,diffmat2{1}(xval,4));
%     end
%     set(gca,'xticklabels',[]);
%     xlabel('ROI'); ylabel('# voxels');
%     title(['mHRF - cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none');
%     set(gca,'xtick',1:length(diffmat2{1}),'xticklabels',roilabels);
%     xtickangle(45)
%     %legend(roilabels,'interpreter','none','Location','NorthEast');
% end
% if SaveFigs; saveas(f3,fullfile(pngfld, 'HRF_Comparison.png')); end
% if SaveFigs; saveas(f3,fullfile(svgfld, 'HRF_Comparison.svg')); end
% if CloseFigs; close(f3); end

f3a=figure;
set(f3a,'Position',[100 100 1300 1200]);
set(f3a,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
s_R2 = T(modidx.linhrf_cv1_mhrf).mod.R2>RTHRES;

f3b=figure;
set(f3b,'Position',[100 100 1300 1200]);
set(f3b,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
s_R2 = T(modidx.linhrf_cv1_mhrf).mod.R2>RTHRES;

for mm = 1:size(MRI_MODELS,1)
    
    figure(f3a)
    subplot(1,4,mm); hold on;
    XY=[];
    for r=1:length(roi)
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                ck_GetROIidx(roilabels(r),rois) );
%         scatter(T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS),...
%             T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS),100,'Marker','.');
        XY=[XY; T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS) ...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS)];
    end
    binscatter(XY(:,1),XY(:,2),100); colorbar; 
    set(gca,'ColorScale','log');
    colormap(inferno)
    caxis([1 1e4]) 
    plot([-2 100],[-2 100],'k','Linewidth',1);   
    title(['mHRF vs cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
    xlabel('Monkey HRF R2'); ylabel('Canonical HRF R2');
    set(gca, 'Box','off', 'xlim', [-2 100], 'ylim',[-2 100]);

    diffmat2{1}=[]; dH=[];
    for r=1:length(roi)
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
                ck_GetROIidx(roilabels(r),rois) );
        [n,x] = hist(...
            T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS)-...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS),100);
        f = n./sum(SSS);
    
        m = mean(T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS)-...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS));
        sd = std(T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS)-...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS));
        nvox = sum(SSS);
        se = sd ./ sqrt(nvox);
        diffmat2{1} = [diffmat2{1}; m sd se nvox];
        dH=[dH;...
            T(modidx.(MRI_MODELS{mm,1})).mod.R2(SSS) - ...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(SSS) ];
    end
    
    figure(f3b)
    subplot(4,1,mm); hold on
    for xval=1:length(diffmat2{1})
        bar(xval,diffmat2{1}(xval,1));
    end
    for xval=1:length(diffmat2{1})
        errorbar(xval,diffmat2{1}(xval,1),diffmat2{1}(xval,3),...
        'k-','Linestyle','none')
    end
    % mean (weighted mean taking nvox into account)
    mAll = sum((diffmat2{1}(:,1).*diffmat2{1}(:,4)))./...
        sum(diffmat2{1}(:,4)); 
    mAll = mean(dH); 
    stdAll = std(dH);
    
    
    % Wilcoxon 1-tailed < 1
    [p,h,stats] = signrank(dH,0);
    fprintf(['ALL VOXELS - Wilcoxon dHRF (' MRI_MODELS{mm,1} ') < 1 : z = ' ...
        num2str(stats.zval) ', p = ' num2str(p) '\n']);
    
    
    
    text(0.5, -1.2,[num2str(mAll) ' +/- ' num2str(stdAll)])
    set(gca,'xticklabels',[],'ylim',[-1.6 1.6],'TickDir','out');
    %xlabel('ROI'); 
    ylabel('Diff R2');
    title(['mHRF - cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
    %legend(roilabels,'interpreter','none','Location','NorthEast');
    set(gca,'xtick',1:length(diffmat2{1}),'xticklabels',roilabels);
    xtickangle(45)    
end
if SaveFigs; saveas(f3b,fullfile(pngfld, 'HRF_Comparison_binned.png')); end
if SaveFigs; saveas(f3b,fullfile(svgfld, 'HRF_Comparison_binned.svg')); end
if CloseFigs; close(f3b); end

%% MRI rf size depending on HRF ===========================================
R2th=5;

f4=figure;
set(f4,'Position',[100 100 1800 1200]);
set(f4,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));

s_R2 = T(modidx.csshrf_cv1_mhrf).mod.R2 > R2th;

for mm = 1:size(MRI_MODELS,1)
    xy_sz{mm}=[]; xy_ecc{mm}=[];
    
    sp=subplot(4,3,(mm-1)*3 +1);hold on;
    plot([0 100],[0 100],'k','LineWidth',2);
    for r=1:length(roi)
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
            ck_GetROIidx(roilabels(r),rois) );
        xy_sz{mm} = [xy_sz{mm};...
            T(modidx.(MRI_MODELS{mm,1})).mod.rfs(SSS) ...
            T(modidx.(MRI_MODELS{mm,2})).mod.rfs(SSS)];
        xy_ecc{mm} = [xy_ecc{mm};...
            T(modidx.(MRI_MODELS{mm,1})).mod.ecc(SSS) ...
            T(modidx.(MRI_MODELS{mm,2})).mod.ecc(SSS)];   
        scatter(T(modidx.(MRI_MODELS{mm,1})).mod.rfs(SSS),...
            T(modidx.(MRI_MODELS{mm,2})).mod.rfs(SSS),100,'Marker','.');
    end
    title(['mHRF vs cHRF (' MMS{mm,1} ')'],'interpreter','none'); 
    xlabel('Monkey HRF sz');ylabel('Canonical HRF sz');
    set(gca, 'Box','off', 'xlim', [0 15], 'ylim',[0 15]);
    text('Parent',sp,'Position',[1 13], ...
        'String',['R2th: ' num2str(R2th)],...
        'Fontsize',12, 'Fontweight','bold')

    xy_sz{mm}( isinf(xy_sz{mm})  ) = nan;
    xy_ecc{mm}( isinf(xy_ecc{mm})  ) = nan;

    % Wilcoxon
    fprintf(['--' MMS{mm,1} '--\n'])
    [p,h,stats] = signrank(xy_sz{mm}(:,1),xy_sz{mm}(:,2),'method','approximate');
    fprintf(['ALL VOXELS R2>TH - Wilcoxon HRF sz: z = ' ...
        num2str(stats.zval) ', p = ' num2str(p) '\n']);
    fprintf(['Mean:' num2str(nanmean(diff(xy_sz{mm},1,2))) ', std ' ...
        num2str(nanstd(diff(xy_sz{mm},1,2))) '\n'])
    
    [p,h,stats] = signrank(xy_ecc{mm}(:,1),xy_ecc{mm}(:,2),'method','approximate');
    fprintf(['ALL VOXELS R2>TH - Wilcoxon HRF ecc: z = ' ...
        num2str(stats.zval) ', p = ' num2str(p) '\n']);
    fprintf(['Mean:' num2str(nanmean(diff(xy_ecc{mm},1,2))) ', std ' ...
        num2str(nanstd(diff(xy_ecc{mm},1,2))) '\n'])
    
    diffmat2{1}=[];
    for r=1:length(roi)
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{rowmod,1})).mod.ROI,...
            ck_GetROIidx(roilabels(r),rois) );
        [n,x] = hist(...
            T(modidx.(MRI_MODELS{mm,1})).mod.rfs(SSS) - ...
            T(modidx.(MRI_MODELS{mm,2})).mod.rfs(SSS), 100);
        f = n./sum(SSS);
        
        dsz = T(modidx.(MRI_MODELS{mm,1})).mod.rfs - ...
            T(modidx.(MRI_MODELS{mm,2})).mod.rfs;
        dsz = dsz(SSS);
        dsz = dsz(isfinite(dsz));
        
        m = mean(dsz);
        sd = std(dsz);
        se = sd ./ sqrt(length(dsz));
        diffmat2{1} = [diffmat2{1}; m sd se size(dsz,1)];
    end
    
    
    subplot(4,3,(mm-1)*3 +2); hold on
    for xval=1:size(diffmat2{1},1)
        bar(xval,diffmat2{1}(xval,1));
    end
    for xval=1:size(diffmat2{1},1)
        errorbar(xval,diffmat2{1}(xval,1),diffmat2{1}(xval,3),...
            'k-','Linestyle','none')
    end
    set(gca,'xticklabels',[],'ylim',[-2 2.5]);
    xlabel('ROI'); ylabel('Diff pRF size');
    title(['mHRF - cHRF (' MMS{mm,1} ')'],'interpreter','none'); 
    %legend(roilabels,'interpreter','none','Location','NorthEast');
    set(gca,'xtick',1:length(diffmat2{1}),'xticklabels',roilabels);
    xtickangle(45)

    subplot(4,3,(mm-1)*3 +3); hold on
    for xval=1:size(diffmat2{1},1)
        bar(xval,diffmat2{1}(xval,4));
    end

    set(gca,'xticklabels',[]);
    xlabel('ROI'); ylabel('Number of voxels');
    title(['mHRF - cHRF (' MMS{mm,1} ')'],'interpreter','none'); 
    set(gca,'xtick',1:length(diffmat2{1}),'xticklabels',roilabels);
    xtickangle(45)
end
if SaveFigs; saveas(f4,fullfile(figfld, 'HRF_Comparison_pRF-Size.png')); end
if CloseFigs; close(f4); end

%% MRI ECC vs PRF Size ====================================================
Rth=5;

f5=figure;
set(f5,'Position',[100 100 1000 1200]);
set(f5,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));

for m=1:length(MRI_MODELS)
    s_R2 = T(modidx.(MRI_MODELS{m,1})).mod.R2 > Rth;
    EccBin = 0.5:1:16.5;
    
    subplot(4,2,(m-1)*2 +1);hold on;
    for r=1:length(roi)
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{m,1})).mod.ROI,...
            ck_GetROIidx(roilabels(r),rois) );
        
        ES{r}=[];
        scatter(T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS),...
            T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS),100,'Marker','.');        
        
        for b=1:length(EccBin)
            bb=[EccBin(b)-0.5 EccBin(b)+0.5];
            PSZ=T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS);
            ECC=T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS);
            ES{r}=[ES{r}; EccBin(b) median(PSZ(ECC>=bb(1) & ECC<=bb(2)))];
        end
    end
    title(['Ecc vs pRF size [' MMS{m} ', R>' num2str(Rth) ']'],...
        'interpreter','none');
    xlabel('Eccentricity');ylabel('pRF size');
    set(gca, 'Box','off', 'xlim', [0 10], 'ylim',[0 10]);
    
    subplot(4,2,(m-1)*2 +2);hold on;
    for r=1:length(roi)
        h=plot(ES{r}(:,1),ES{r}(:,2),'o');
        set(h,'MarkerSize',6,'markerfacecolor', get(h, 'color'));
    end
    title(['Ecc vs pRF size [' MMS{m} ', R>' num2str(Rth) ']'],...
        'interpreter','none'); xlabel('Eccentricity');ylabel('pRF size');
    set(gca, 'Box','off', 'xlim', [0 10], 'ylim',[0 10]);
    %legend(roilabels,'interpreter','none','Location','NorthWest');
end
if SaveFigs
    saveas(f5,fullfile(pngfld, 'MRI_Ecc_vs_Size.png'));
    saveas(f5,fullfile(svgfld, 'MRI_Ecc_vs_Size.svg'));
end
if CloseFigs; close(f5); end

%% More detail for CSS mHRF (for paper) -----------------------------------
Rth=5;
maxEcc = 50;

f5css=figure;
set(f5css,'Position',[100 100 2000 1500]);
set(f5css,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));

m=3; % css
s_R2 = T(modidx.(MRI_MODELS{m,1})).mod.R2 > Rth;

EccBin = 0.5:1:16.5;

nsp=length(roi);
nrc=ceil(sqrt(nsp));

for r=1:length(roi)
    ES{r} = []; ES_m{r} = [];
    ESb{r} = [];
    ESf{r} = [];
    
    SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{m,1})).mod.ROI,...
        ck_GetROIidx(roilabels(r),rois) );
    s_Monkey = strcmp(T(modidx.(MRI_MODELS{m,1})).mod.Monkey,'danny');
    SSS_m1 = s_R2 & ismember( T(modidx.(MRI_MODELS{m,1})).mod.ROI,...
        ck_GetROIidx(roilabels(r),rois)) & s_Monkey;
    s_Monkey = strcmp(T(modidx.(MRI_MODELS{m,1})).mod.Monkey,'eddy');
    SSS_m2 = s_R2 & ismember( T(modidx.(MRI_MODELS{m,1})).mod.ROI,...
        ck_GetROIidx(roilabels(r),rois)) & s_Monkey;
    
    ES{r}=[T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS) ...
        T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS) ];
    ES{r}( ES{r}(:,2) > 50,: ) = []; % remove insanely large pRF's
    
    ES_m{1,r}=[T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS_m1) ...
        T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS_m1) ];
    ES_m{2,r}=[T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS_m2) ...
        T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS_m2) ];
    
    subplot(nrc,nrc,r); hold on;
    

%     scatter(ES{r}(:,1),ES{r}(:,2),10,'Marker','o',...
%             'MarkerFaceColor','k','MarkerFaceAlpha',.2,...
%             'MarkerEdgeColor','none');
    scatter(ES_m{1,r}(:,1),ES_m{1,r}(:,2),10,'Marker','o',...
            'MarkerFaceColor','b','MarkerFaceAlpha',.2,...
            'MarkerEdgeColor','none');
    scatter(ES_m{2,r}(:,1),ES_m{2,r}(:,2),10,'Marker','o',...
            'MarkerFaceColor','k','MarkerFaceAlpha',.2,...
            'MarkerEdgeColor','none');
        
    % make table 
    tbl = table(ES{r}(:,1),ES{r}(:,2),'VariableNames',{'ecc','sz'});
    if ~isempty(tbl)
        % remove inf's
        tbl(isinf(tbl.sz),:)=[];
        % fit
        lm{r} = fitlm(tbl,'sz ~ 1 + ecc');
        % CI
        lmCI{r}=coefCI(lm{r},0.05);
        if ~isnan(lm{r}.Coefficients.pValue(1))
            % Prediction
            x=[0;20];
            [ypred,yci] = predict(lm{r},x);
            yci2=yci(:,2)-ypred;
            boundedline(x,ypred',yci2,'r','alpha')
            %fprintf([roilabels{r} ' n = ' num2str(size(tbl,1)) '\n'])
            %fprintf([roilabels{r} ' p = ' num2str(lm{r}.Coefficients.pValue(2)) '\n'])
            fprintf([roilabels{r} ' n = ' num2str(lm{r}.Coefficients.Estimate(2)) '\n'])
        end
    end
    
%     for b=1:length(EccBin)
%         bb=[EccBin(b)-0.5 EccBin(b)+0.5];
%         PSZ=T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS);
%         ECC=T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS);
%         ES{r}=[ES{r}; EccBin(b) median(PSZ(ECC>=bb(1) & ECC<=bb(2)))];
%     end
    title(['Ecc vs pRF size [' roilabels{r} ', R>' num2str(Rth) ']'],...
    'interpreter','none');
    title(roilabels{r});
    xlabel('Eccentricity');ylabel('pRF size');
    set(gca, 'Box','off', 'xlim', [0 20], 'ylim',[0 25]);
end

if SaveFigs
    saveas(f5css,fullfile(pngfld, 'MRI_Ecc_vs_Size.png'));
    saveas(f5css,fullfile(svgfld, 'MRI_Ecc_vs_Size.svg'));
end
if CloseFigs; close(f5css); end

%% Create manuscript figure -----------------------------------------------
subctx_idx = [18:20];
early_ctx_idx = [1:5 7:14];
late_ctx_idx = [23 25:27 31 33:35 38];

f_eccsz = figure;
set(f_eccsz,'Position',[100 100 2000 1500]);

subplot(1,3,1); hold on;
lidx=1;
cm=brewermap(length(subctx_idx),def_cmap);
set(f_eccsz,'DefaultAxesColorOrder',cm);
for idx=subctx_idx
    x=[0;20];
    [ypred,yci] = predict(lm{idx},x);
    yci2=yci(:,2)-ypred;
    plot(x,ypred','LineWidth',3)
    L{lidx}=roilabels{idx}; lidx=lidx+1;
end
for idx=subctx_idx
    x=[0;20];
    [ypred,yci] = predict(lm{idx},x);
    yci2=yci(:,2)-ypred;
    boundedline(x,ypred',yci2,'o','alpha')
end
legend(L)
set(gca,'ylim',[0 30], 'TickDir','out');


subplot(1,3,2); hold on;
lidx=1;
cm=brewermap(length(early_ctx_idx),def_cmap);
set(f_eccsz,'DefaultAxesColorOrder',cm);
for idx=early_ctx_idx
    x=[0;20];
    [ypred,yci] = predict(lm{idx},x);
    yci2=yci(:,2)-ypred;
    plot(x,ypred','LineWidth',3)
    L{lidx}=roilabels{idx}; lidx=lidx+1;
end
for idx=early_ctx_idx
    x=[0;20];
    [ypred,yci] = predict(lm{idx},x);
    yci2=yci(:,2)-ypred;
    boundedline(x,ypred',yci2,'o','alpha')
end
legend(L)
set(gca,'ylim',[0 30], 'TickDir','out');

subplot(1,3,3); hold on;
lidx=1;
cm=brewermap(length(late_ctx_idx),def_cmap);
set(f_eccsz,'DefaultAxesColorOrder',cm);
for idx=late_ctx_idx
    x=[0;20];
    [ypred,yci] = predict(lm{idx},x);
    yci2=yci(:,2)-ypred;
    plot(x,ypred','LineWidth',3)
    L{lidx}=roilabels{idx}; lidx=lidx+1;
end
for idx=late_ctx_idx
    x=[0;20];
    [ypred,yci] = predict(lm{idx},x);
    yci2=yci(:,2)-ypred;
    boundedline(x,ypred',yci2,'o','alpha')
end
legend(L)
set(gca,'ylim',[0 30], 'TickDir','out');

%% ------------------------------------------------------------------------
f5b=figure;
set(f5b,'Position',[100 100 1000 1200]);
set(f5b,'DefaultAxesColorOrder',brewermap(2,def_cmap));

for m=1:length(MRI_MODELS)
    s_R2 = T(modidx.(MRI_MODELS{m,1})).mod.R2 > Rth;
    EccBin = 0.5:1:20.5;
    
    subplot(4,2,(m-1)*2 +1);hold on;
    for r=[1 4]
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{m,1})).mod.ROI,...
            ck_GetROIidx(roilabels(r),rois) );

        ES{r}=[];
        scatter(T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS),...
            T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS),100,'Marker','.');
        for b=1:length(EccBin)
            bb=[EccBin(b)-0.5 EccBin(b)+0.5];
            PSZ=T(modidx.(MRI_MODELS{m,1})).mod.rfs(SSS);
            ECC=T(modidx.(MRI_MODELS{m,1})).mod.ecc(SSS);
            ES{r}=[ES{r}; EccBin(b) median(PSZ(ECC>=bb(1) & ECC<=bb(2)))];
        end
    end
    title(['Ecc vs pRF size [' MMS{m} ', R>' num2str(Rth) ']'],...
        'interpreter','none');
    xlabel('Eccentricity');ylabel('pRF size');
    set(gca, 'Box','off', 'xlim', [0 10], 'ylim',[0 10]);
    
    subplot(4,2,(m-1)*2 +2);hold on;
    for r=[1 4]
        h=plot(ES{r}(:,1),ES{r}(:,2),'o');
        set(h,'MarkerSize',6,'markerfacecolor', get(h, 'color'));
    end
    title(['Ecc vs pRF size [' MMS{m} ', R>' num2str(Rth) ']'],...
        'interpreter','none'); xlabel('Eccentricity');ylabel('pRF size');
    set(gca, 'Box','off', 'xlim', [0 10], 'ylim',[0 5]);
    %legend(roilabels,'interpreter','none','Location','NorthWest');
end
if SaveFigs; saveas(f5b,fullfile(figfld, 'MRI_Ecc_vs_Size_V1V4.png')); end
if CloseFigs; close(f5b); end






%% ========================================================================
%  EPHYS PRF ANALYSIS -----------------------------------------------------
%  ========================================================================

%% EPHYS ECC vs PRF Size ==================================================
Rth=25; 
SNRth=2;

m=unique(tMUA.Model);
R2=[];
for i=1:length(m)
    R2 = [R2 tMUA.R2(strcmp(tMUA.Model,m{i}))];
end

subs = tMUA.Monkey;
for i=1:length(subs)
    if strcmp(subs{i},'lick')
        subs{i} = 1;
    elseif strcmp(subs{i},'aston')
        subs{i} = 2;
    end
end
subs=cell2mat(subs); names={'lick','aston'};

mm=unique(tMUA.Model);
for m=1:length(mm)
    % collect
    if ~strcmp(mm{m},'classicRF')
        eccsz = [tMUA.R2(strcmp(tMUA.Model,mm{m})) ...
            tMUA.ecc(strcmp(tMUA.Model,mm{m})) ...
            tMUA.rfs(strcmp(tMUA.Model,mm{m})) ...
            tMUA.Area(strcmp(tMUA.Model,mm{m})) ...
            tMUA.Array(strcmp(tMUA.Model,mm{m})) ...
            subs(strcmp(tMUA.Model,mm{m})) ...
            tMUA.gain(strcmp(tMUA.Model,mm{m})) ...
            ];
    else 
        eccsz = [tMUA.SNR(strcmp(tMUA.Model,mm{m})) ...
            tMUA.ecc(strcmp(tMUA.Model,mm{m})) ...
            tMUA.rfs(strcmp(tMUA.Model,mm{m})) ...
            tMUA.Area(strcmp(tMUA.Model,mm{m})) ...
            tMUA.Array(strcmp(tMUA.Model,mm{m})) ...
            subs(strcmp(tMUA.Model,mm{m})) ...
            ];
    end
    
    f_eccsz_arr = figure;
    set(f_eccsz_arr,'Position',[100 100 1800 600]);
    for ss=1:2
        a=[1 4];
        for aa=1:length(a)
            leg={};
            if strcmp(mm{m},'classicRF')
                sel = eccsz(:,1)>SNRth & eccsz(:,4)==a(aa) & eccsz(:,6)==ss;
            else
                sel = eccsz(:,1)>Rth & eccsz(:,4)==a(aa) & eccsz(:,6)==ss;
            end
            es1=eccsz(sel,:);
            subplot(2,2,(ss-1)*2 + aa); 
            hold on;
            arr=unique(es1(:,5));
            for aaa = 1:length(arr)
                if arr(aaa)~=99
                    sel2=es1(:,5)==arr(aaa);
                    scatter(es1(sel2,2),es1(sel2,3));
                    leg=[leg {num2str(arr(aaa))}];
                end
            end
            if a(aa)==1
                set(gca,'xlim',[0 8],'ylim',[0 3]);
            else
                set(gca,'xlim',[0 8],'ylim',[0 5]);
            end
            legend(leg,'Location','NorthEastOutside');
            title([names{ss} ' V' num2str(a(aa))]);
            xlabel('Ecc (dva)'); ylabel('RF-size (sd, dva)')
        end
    end
    sgtitle(['MODEL: ' mm{m} ', R2th: ' num2str(Rth)],...
        'interpreter','none');
    if SaveFigs
        saveas(f_eccsz_arr,fullfile(figfld, ...
            ['EPHYS_MUA_Ecc_vs_Size_' mm{m} '.png']));
    end
    if CloseFigs; close(f_eccsz_arr); end
    
    f_eccsz_v1 = figure;
    set(f_eccsz_v1,'Position',[100 100 1600 400]);
    for ss=1:2
        leg={};
        if strcmp(mm{m},'classicRF')
            sel = eccsz(:,1)>SNRth & eccsz(:,4)==1 & eccsz(:,6)==ss;
        else
            sel = eccsz(:,1)>Rth & eccsz(:,4)==1 & eccsz(:,6)==ss;
        end
        es1=eccsz(sel,:);
        subplot(1,2,ss); hold on;
        arr=unique(es1(:,5));
        for aaa = 1:length(arr)
            sel2=es1(:,5)==arr(aaa);
            scat=scatter(es1(sel2,2),es1(sel2,3));
            %if (ss==1 && (arr(aaa)==11 || arr(aaa)==13)) || ...
            %        (ss==2 && (arr(aaa)==10 || arr(aaa)==12))
            if (ss==1 && arr(aaa)==11) || ...
                    (ss==2 && arr(aaa)==10)
                set(scat,'MarkerFaceColor','r','MarkerEdgeColor','r');
            elseif (ss==1 && arr(aaa)==13) || ...
                    (ss==2 && arr(aaa)==12)
                set(scat,'MarkerFaceColor','b','MarkerEdgeColor','b');
                
            else
                set(scat,'MarkerFaceColor','k','MarkerEdgeColor','k');
            end
            leg=[leg {num2str(arr(aaa))}];
        end
        title([names{ss} ' V1']);
        xlabel('Ecc (dva)'); ylabel('RF-size (sd, dva)')
        set(gca,'xlim',[0 8],'ylim',[0 3]);
        %legend(leg,'Location','NorthEastOutside');
    end
    sgtitle(['MODEL: ' mm{m} ', R2th: ' num2str(Rth)],...
        'interpreter','none');
    if SaveFigs
        saveas(f_eccsz_v1,fullfile(figfld, ...
            ['EPHYS_MUA_Ecc_vs_Size_' mm{m} '.png']));
    end
    if CloseFigs; close(f_eccsz_v1); end
    
    % ====================================
    % there's something odd about arrays:
    % 11 [to a lesser extent 13] in lick
    % 10 [to a lesser extent 12] in aston
    % ====================================

%     % check their SNR (from classic rf fits)
%     % these seem pretty normal
%     figure;
%     subplot(1,2,1);
%     ss=eccsz(:,6)==1;
%     scatter(eccsz(ss,5),eccsz(ss,1))
%     subplot(1,2,2);
%     ss=eccsz(:,6)==2;
%     scatter(eccsz(ss,5),eccsz(ss,7))
%     
%     % check their gains
%     % also normal
%     figure;
%     subplot(1,2,1);
%     ss=eccsz(:,6)==1;
%     scatter(eccsz(ss,5),eccsz(ss,7))
%     subplot(1,2,2);
%     ss=eccsz(:,6)==2;
%     scatter(eccsz(ss,5),eccsz(ss,7))
   
end
%close all;


% focus on CSS --------
% exclude the problematic arrays

%% EPHYS VFC MUA ==========================================================

% scatter MUA =============================================================
R2TH = 50;
s=tMUA.R2>R2TH;
model='css_ephys_cv1';

fm=figure; set(fm,'Position',[100 100 1800 1400]);
[cmap,~,~] = brewermap(length(unique(tMUA.Array)),'spectral');
set(fm,'DefaultAxesColorOrder',cmap);
msz = 5;

% Lick V1 ----
m=strcmp(tMUA.Monkey,'lick') & strcmp(tMUA.Model,model);
v=tMUA.Area==1;

subplot(2,2,1);hold on;
plot([-10 20],[0 0],'k'); plot([0 0],[-30 10],'k');
lt={'meridian','meridian'};
nelec=0;
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    currcol = cmap(r,:);
    if sum(s & m & a & v) > 0
        nelec=nelec+sum(s & m & a & v);
        p=plot(tMUA.X(s & m & a & v),...
            tMUA.Y(s & m & a & v),'o','Color',currcol,...
            'LineStyle','none','LineWidth',1,...
            'MarkerSize',msz,'MarkerFaceColor',currcol);
        lt=[lt num2str(r)];
    end
end
set(gca, 'Box','off', 'xlim', [-1 8], 'ylim',[-6 1]);
title(['Lick V1 MUA (n = ' num2str(nelec) ')'])
l=legend(lt); set(l,'Location','NorthEastOutside'); clear lt;

% Lick V4 ----
v=tMUA.Area==4;

subplot(2,2,3);hold on;
plot([-10 20],[0 0],'k'); plot([0 0],[-30 10],'k');
lt={'meridian','meridian'};
nelec=0;
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    currcol = cmap(r,:);
    if sum(s & m & a & v) > 0
        nelec=nelec+sum(s & m & a & v);
        p=plot(tMUA.X(s & m & a & v),...
            tMUA.Y(s & m & a & v),'o','Color',currcol,...
            'LineStyle','none','LineWidth',1,...
            'MarkerSize',msz,'MarkerFaceColor',currcol);
        lt=[lt num2str(r)];
    end
end
set(gca, 'Box','off', 'xlim', [-1 8], 'ylim',[-6 1]);
% set(gca, 'Box','off', 'xlim', [-2 25], 'ylim',[-25 2]);
title(['Lick V4 MUA (n = ' num2str(nelec) ')'])
l=legend(lt); set(l,'Location','NorthEastOutside'); clear lt;

% Aston V1 ----
m=strcmp(tMUA.Monkey,'aston') & strcmp(tMUA.Model,model);
v=tMUA.Area==1;

subplot(2,2,2);hold on;
plot([-10 20],[0 0],'k'); plot([0 0],[-30 10],'k');
lt={'meridian','meridian'};
nelec=0;
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    currcol = cmap(r,:);
    if sum(s & m & a & v) > 0
        nelec=nelec+sum(s & m & a & v);
        p=plot(tMUA.X(s & m & a & v),...
            tMUA.Y(s & m & a & v),'o','Color',currcol,...
            'LineStyle','none','LineWidth',1,...
            'MarkerSize',msz,'MarkerFaceColor',currcol);
        lt=[lt num2str(r)];
    end
end
set(gca, 'Box','off', 'xlim', [-1 8], 'ylim',[-6 1]);
title(['Aston V1 MUA (n = ' num2str(nelec) ')'])
l=legend(lt); set(l,'Location','NorthEastOutside'); clear lt;

% Aston V4 ----
v=tMUA.Area==4;

subplot(2,2,4);hold on;
plot([-10 30],[0 0],'k'); plot([0 0],[-30 10],'k');
lt={'meridian','meridian'};
nelec=0;
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    currcol = cmap(r,:);
    if sum(s & m & a & v) > 0
        nelec=nelec+sum(s & m & a & v);
        p=plot(tMUA.X(s & m & a & v),...
            tMUA.Y(s & m & a & v),'o','Color',currcol,...
            'LineStyle','none','LineWidth',1,...
            'MarkerSize',msz,'MarkerFaceColor',currcol);
        lt=[lt num2str(r)];
    end
end
set(gca, 'Box','off', 'xlim', [-2 25], 'ylim',[-25 2]);
title(['Aston V4 MUA (n = ' num2str(nelec) ')'])
l=legend(lt); set(l,'Location','NorthEastOutside'); clear lt;

% Heatmap MUA =============================================================
fhm=figure; set(fhm,'Position',[100 100 1800 1000]);
settings.PixPerDeg = 29.5032;
settings.meshsize = 2000;   
colormap(inferno) 

% Lick V1 ----
m=strcmp(tMUA.Monkey,'lick') & strcmp(tMUA.Model,model);
v=tMUA.Area==1;

subplot(2,4,1);hold on;
allprf=[];
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    if sum(s & m & a & v) > 0
        allprf = [allprf; ...
            tMUA.X(s & m & a & v) ...
            tMUA.Y(s & m & a & v) ...
            tMUA.rfs(s & m & a & v)];
    end
end

res = ck_2dPRF_ephys(allprf(:,1),allprf(:,2),allprf(:,3));
img=flipud(res.img); res.ymesh2 = fliplr(res.ymesh);
% zoom in on plot
xrange = [-3 8]; yrange = [-6 3];
xrange_idx = [...
    find(res.xmesh >= xrange(1),1,'first') find(res.xmesh <= xrange(2),1,'last')];
yrange_idx = [...
    find(res.ymesh2 >= yrange(1),1,'first') find(res.ymesh2 <= yrange(2),1,'last')];
xrr=res.xmesh(xrange_idx); yrr=res.ymesh2(yrange_idx);
% plot
sumimg=sum(img,3);
imagesc(sumimg);
set(gca,'xlim',xrange_idx,'ylim',...
    yrange_idx,'Color','k','xtick',[],'ytick',[])
colorbar;
subplot(2,4,5); hold on;
plot(1.2*res.xr,[0 0],'w'); plot([0 0],1.2*res.yr,'w');
set(gca,'xlim',xrr,'ylim',yrr,'Color','k')
colorbar; clear('res','img','sumimg');
title('Lick V1 MUA')

% Lick V4 ----
v=tMUA.Area==4;

subplot(2,4,2);hold on;
allprf=[];
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    if sum(s & m & a & v) > 0
        allprf = [allprf; ...
            tMUA.X(s & m & a & v) ...
            tMUA.Y(s & m & a & v) ...
            tMUA.rfs(s & m & a & v)];
    end
end

res = ck_2dPRF_ephys(allprf(:,1),allprf(:,2),allprf(:,3));
img=flipud(res.img); res.ymesh2 = fliplr(res.ymesh);
% zoom in on plot
xrange = [-3 8]; yrange = [-6 3];
xrange_idx = [...
    find(res.xmesh >= xrange(1),1,'first') find(res.xmesh <= xrange(2),1,'last')];
yrange_idx = [...
    find(res.ymesh2 >= yrange(1),1,'first') find(res.ymesh2 <= yrange(2),1,'last')];
xrr=res.xmesh(xrange_idx); yrr=res.ymesh2(yrange_idx);
% plot
sumimg=sum(img,3);
imagesc(sumimg);
set(gca,'xlim',xrange_idx,'ylim',...
    yrange_idx,'Color','k','xtick',[],'ytick',[])
colorbar;
subplot(2,4,6); hold on;
plot(1.2*res.xr,[0 0],'w'); plot([0 0],1.2*res.yr,'w');
set(gca,'xlim',xrr,'ylim',yrr,'Color','k')
colorbar; clear('res','img','sumimg');
title('Lick V4 MUA')

% Aston V1 ----
m=strcmp(tMUA.Monkey,'aston') & strcmp(tMUA.Model,model);
v=tMUA.Area==1;

subplot(2,4,3);hold on;
allprf=[];
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    if sum(s & m & a & v) > 0
        allprf = [allprf; ...
            tMUA.X(s & m & a & v) ...
            tMUA.Y(s & m & a & v) ...
            tMUA.rfs(s & m & a & v)];
    end
end

res = ck_2dPRF_ephys(allprf(:,1),allprf(:,2),allprf(:,3));
img=flipud(res.img); res.ymesh2 = fliplr(res.ymesh);
% zoom in on plot
xrange = [-3 8]; yrange = [-6 3];
xrange_idx = [...
    find(res.xmesh >= xrange(1),1,'first') find(res.xmesh <= xrange(2),1,'last')];
yrange_idx = [...
    find(res.ymesh2 >= yrange(1),1,'first') find(res.ymesh2 <= yrange(2),1,'last')];
xrr=res.xmesh(xrange_idx); yrr=res.ymesh2(yrange_idx);
% plot
sumimg=sum(img,3);
imagesc(sumimg);
set(gca,'xlim',xrange_idx,'ylim',...
    yrange_idx,'Color','k','xtick',[],'ytick',[])
colorbar;
subplot(2,4,7); hold on;
plot(1.2*res.xr,[0 0],'w'); plot([0 0],1.2*res.yr,'w');
set(gca,'xlim',xrr,'ylim',yrr,'Color','k')
colorbar; clear('res','img','sumimg');
title('Aston V1 MUA')

% Aston V4 ----
v=tMUA.Area==4;

subplot(2,4,4);hold on;
allprf=[];
for r=unique(tMUA.Array)'
    a=tMUA.Array==r;
    if sum(s & m & a & v) > 0
        allprf = [allprf; ...
            tMUA.X(s & m & a & v) ...
            tMUA.Y(s & m & a & v) ...
            tMUA.rfs(s & m & a & v)];
    end
end

res = ck_2dPRF_ephys(allprf(:,1),allprf(:,2),allprf(:,3));
img=flipud(res.img); res.ymesh2 = fliplr(res.ymesh);
% zoom in on plot
xrange = [-5 25]; yrange = [-25 5];
xrange_idx = [...
    find(res.xmesh >= xrange(1),1,'first') find(res.xmesh <= xrange(2),1,'last')];
yrange_idx = [...
    find(res.ymesh2 >= yrange(1),1,'first') find(res.ymesh2 <= yrange(2),1,'last')];
xrr=res.xmesh(xrange_idx); yrr=res.ymesh2(yrange_idx);
% plot
sumimg=sum(img,3);
imagesc(sumimg);
set(gca,'xlim',xrange_idx,'ylim',...
    yrange_idx,'Color','k','xtick',[],'ytick',[])
colorbar;
subplot(2,4,8); hold on;
plot(1.2*res.xr,[0 0],'w'); plot([0 0],1.2*res.yr,'w');
set(gca,'xlim',xrr,'ylim',yrr,'Color','k')
colorbar; clear('res','img','sumimg');
title('Aston V4 MUA')

%% Ephys VFC LFP GAMMA ====================================================
% LFP Low Gamma
s=tLFP.R2>50;

model='css_ephys_cv1';
b = 'lGamma';
fm=figure;
set(fm,'Position',[100 100 800 800]);

subplot(2,2,1);hold on;
szm=tLFP.rfs<50;
m=strcmp(tLFP.Monkey,'lick') & ...
    strcmp(tLFP.Model, model) & ...
    strcmp(tLFP.SigType, b);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP.Area==1;
for r=unique(tLFP.Array)'
    a=tLFP.Array==r;
    scatter(tLFP.X(s & m & a & v),...
        tLFP.Y(s & m & a & v),'Marker','*' )
end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Lick V1 LFP (lGAM)')

subplot(2,2,3);hold on;
szm=tLFP.rfs<50;
m=strcmp(tLFP.Monkey,'lick') & strcmp(tLFP.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP.Area==4;
for r=unique(tLFP.Array)'
    a=tLFP.Array==r;
    scatter(tLFP.X(s & m & a & v),...
        tLFP.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Lick V4 LFP (lGAM)')

subplot(2,2,2);hold on;
m=strcmp(tLFP.Monkey,'aston') & strcmp(tLFP.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP.Area==1;
for r=unique(tLFP.Array)'
    a=tLFP.Array==r;
    scatter(tLFP.X(s & m & a & v),...
        tLFP.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Aston V1 LFP (lGAM)')

subplot(2,2,4);hold on;
m=strcmp(tLFP.Monkey,'aston') & strcmp(tLFP.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP.Area==4;
for r=unique(tLFP.Array)'
    a=tLFP.Array==r;
    scatter(tLFP.X(s & m & a & v),...
        tLFP.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Aston V4 LFP (lGAM)')

if SaveFigs
    saveas(fm,fullfile(figfld, ['EPHYS_VFC_LFP_' model '.png']));
end
if CloseFigs; close(fm); end

%% Ephys location difference & size difference  ---------------------------
rth=25; snrth=3;

ephys_MMS = MMS(:,1);
clear C R2m SZ
for m=1:length(ephys_MOD)
 
    C{m}=[];R2m{m}=[];SZ{m}=[];
    
    model=ephys_MOD{m};
    s = strcmp(tMUA.Model,model);
    C{m}=[C{m} tMUA.R2(s) tMUA.X(s) tMUA.Y(s) tMUA.rfs(s)];
    R2m{m}=[R2m{m} tMUA.R2(s)];
    SZ{m}=[SZ{m} tMUA.R2(s) tMUA.rfs(s)];
    
    PRF_EST(m,1).sig = 'MUA';
    PRF_EST(m,1).R2 = tMUA.R2(s);
    PRF_EST(m,1).X = tMUA.X(s);
    PRF_EST(m,1).Y = tMUA.Y(s);
    PRF_EST(m,1).S = tMUA.rfs(s);
    PRF_EST(m,1).A = tMUA.Area(s);
    PRF_EST(m,1).G = tMUA.gain(s);
    PRF_EST(m,1).M = tMUA.Monkey(s);
    
    s = strcmp(tMUA.Model,'classicRF');
    C{m}=[C{m} tMUA.X(s)./668.745 tMUA.Y(s)./668.745 tMUA.rfs(s)./2];
    
    PRF_EST(m,2).sig = 'MUACLASSIC';
    PRF_EST(m,2).R2 = tMUA.SNR(s);
    PRF_EST(m,2).X = tMUA.X(s);
    PRF_EST(m,2).Y = tMUA.Y(s);
    PRF_EST(m,2).S = tMUA.rfs(s)./2;
    PRF_EST(m,2).A = tMUA.Area(s);
    PRF_EST(m,2).G = tMUA.gain(s);
    PRF_EST(m,2).M = tMUA.Monkey(s);
    
    s = strcmp(tLFP.Model,model);
    sig=unique(tLFP.SigType);
    lfp_order = [3 1 2 5 4];
    cnt=1;
    for i=lfp_order
        b=strcmp(tLFP.SigType,sig{i});
        C{m}=[C{m} tLFP.R2(s & b) tLFP.X(s & b) tLFP.Y(s & b) tLFP.rfs(s & b)];
        R2m{m}=[R2m{m} tLFP.R2(s & b)];
        SZ{m}=[SZ{m} tLFP.R2(s & b) tLFP.rfs(s & b)];
        
        PRF_EST(m,2+cnt).sig = sig{i};
        PRF_EST(m,2+cnt).R2 = tLFP.R2(s & b);
        PRF_EST(m,2+cnt).X = tLFP.X(s & b);
        PRF_EST(m,2+cnt).Y = tLFP.Y(s & b);
        PRF_EST(m,2+cnt).S =  tLFP.rfs(s & b);
        PRF_EST(m,2+cnt).A =  tLFP.Area(s & b);
        PRF_EST(m,2+cnt).G =  tLFP.gain(s & b);
        
        PRF_EST(m,2+cnt).M = tLFP.Monkey(s & b);
        
        cnt=cnt+1;
    end
    
    s= sum(R2m{m}>rth,2)==size(R2m{m},2);
    distRF{m} = [...
        sqrt(((C{m}(s,2)-C{m}(s,5)).^2) + ((C{m}(s,3)-C{m}(s,6)).^2)) ...
        sqrt(((C{m}(s,2)-C{m}(s,9)).^2) + ((C{m}(s,3)-C{m}(s,10)).^2)) ...
        sqrt(((C{m}(s,2)-C{m}(s,13)).^2) + ((C{m}(s,3)-C{m}(s,14)).^2)) ...
        sqrt(((C{m}(s,2)-C{m}(s,17)).^2) + ((C{m}(s,3)-C{m}(s,18)).^2)) ...
        sqrt(((C{m}(s,2)-C{m}(s,21)).^2) + ((C{m}(s,3)-C{m}(s,22)).^2)) ];
    
    distSZ{m} = [...
        C{m}(s,4)-C{m}(s,7) ...
        C{m}(s,4)-C{m}(s,11) ...
        C{m}(s,4)-C{m}(s,15) ...
        C{m}(s,4)-C{m}(s,19) ...
        C{m}(s,4)-C{m}(s,23) ];
       
    normSz{m} = [...
        C{m}(s,7)./C{m}(s,4) ...
        C{m}(s,11)./C{m}(s,4) ...
        C{m}(s,15)./C{m}(s,4) ...
        C{m}(s,19)./C{m}(s,4) ...
        C{m}(s,23)./C{m}(s,4) ];
end

figure;
% compare MUA with classic
for m=3 %1:length(ephys_MOD)
    pn=0;
    for area = [1 4]
        pn=pn+1;
        fprintf(['=== AREA V' num2str(area) ' ===\n'])
        % location
        sel = PRF_EST(m,1).A == area & ...
            PRF_EST(m,1).R2 > rth & ...
            PRF_EST(m,2).R2 >= snrth;
        dLOCATION = sqrt(...
            (PRF_EST(m,1).X(sel)-PRF_EST(m,2).X(sel)).^2 + ...
            (PRF_EST(m,1).Y(sel)-PRF_EST(m,2).Y(sel)).^2);
        
        mselect = PRF_EST(m,1).M(sel);
        m3idx = strcmp(mselect,'aston'); m4idx = strcmp(mselect,'lick');
        
        fprintf(['MODEL ' ephys_MOD{m} ', MUA vs CLASSIC distance ----\n'])
        fprintf(['Mean ' num2str(mean(dLOCATION)) ', STD ' num2str(std(dLOCATION)) '\n'])
        fprintf(['Median ' num2str(median(dLOCATION)) ', IQR ' num2str(iqr(dLOCATION)) '\n'])
        %     figure;
        %     subplot(1,2,1)
        %     scatter(PRF_EST(m,1).X(sel),PRF_EST(m,2).X(sel))
        %     subplot(1,2,2)
        %     scatter(PRF_EST(m,1).Y(sel),PRF_EST(m,2).Y(sel))
        % size
        SIZE_MUA = [PRF_EST(m,1).S(sel) PRF_EST(m,2).S(sel)];
        rmINF = isinf(SIZE_MUA(:,1));
        SIZE_MUA(rmINF,:)=[];
        m3idx(rmINF,:)=[];
        m4idx(rmINF,:)=[];
        
        [p,h,stats] = signrank(SIZE_MUA(:,1),SIZE_MUA(:,2));
        fprintf(['MODEL ' ephys_MOD{m} ', MUA vs CLASSIC size ----\n'])
        fprintf(['dSZ: Wilcoxon z = ' ...
            num2str(stats.zval) ', p = ' num2str(p) '\n']);
        fprintf(['Mean (mua-classic) ' num2str(nanmean(SIZE_MUA(:,1)-SIZE_MUA(:,2))) ...
            ', STD ' num2str(nanstd(SIZE_MUA(:,1)-SIZE_MUA(:,2))) '\n'])
        fprintf(['Median ' num2str(median(SIZE_MUA(:,1)-SIZE_MUA(:,2))) ...
            ', IQR ' num2str(iqr(SIZE_MUA(:,1)-SIZE_MUA(:,2))) '\n'])
        subplot(1,2,pn);hold on;
        splim=[3 8];
        plot([0 10],[0 10],'-r');
%         % all
%         scatter(SIZE_MUA(:,1),SIZE_MUA(:,2),50,'Marker','o',...
%             'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.3 .3 .3],...
%             'MarkerFaceAlpha',0.5);
        % split by monkey
        scatter(SIZE_MUA(m3idx,1),SIZE_MUA(m3idx,2),50,'Marker','o',...
            'MarkerEdgeColor','none',...
            'MarkerFaceColor',[.0 .0 .5],'MarkerFaceAlpha',0.3);
        polyfit(SIZE_MUA(m3idx,1),SIZE_MUA(m3idx,2),1)
        
        
        
        scatter(SIZE_MUA(m4idx,1),SIZE_MUA(m4idx,2),50,'Marker','o',...
            'MarkerEdgeColor','none',...
            'MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',0.3);
        
        plot([0 10],[0 10],'-r');
        set(gca, 'xlim',[0 splim(pn)],'ylim',[0 splim(pn)]);
        title(['Area V' num2str(area)])
        xlabel('MUA pRF size'); ylabel('Classic RF size')
    end
end

%% MUA model comparison ===================================================
m=unique(tMUA.Model);
R2=[];
RTH = 0;

for i=1:length(m)
    R2 = [R2 tMUA.R2(strcmp(tMUA.Model,m{i}))];
end

v1=tMUA.Area(strcmp(tMUA.Model,m{1}))==1;
v4=tMUA.Area(strcmp(tMUA.Model,m{1}))==4;

f6=figure;
msz=15;
set(f6,'Position',[100 100 1200 1200]);
for row=1:4
    for column=1:4
        subplot(4,4,((row-1)*4)+column); hold on;
        plot([0 100],[0 100],'k');
        
        XY = [R2(v1,row+1), R2(v1,column+1)];
        XYs = XY(XY(:,1)>RTH & XY(:,2)>RTH,:);
        scatter(XYs(:,1), XYs(:,2),msz,'Marker','o',...
            'MarkerEdgeColor',[.3 .3 .3],'MarkerFaceColor',[.3 .3 .3],...
            'MarkerFaceAlpha',0.5);
        XY = [R2(v4,row+1), R2(v4,column+1)];
        XYs = XY(XY(:,1)>RTH & XY(:,2)>RTH,:);
        scatter(XYs(:,1), XYs(:,2),msz,'Marker','o',...
            'MarkerEdgeColor',[.3 .8 .3],'MarkerFaceColor',[.3 .8 .3],...
            'MarkerFaceAlpha',0.5);
        set(gca, 'Box','off', 'xlim', [-2 100], 'ylim',[-2 100]);
        xlabel(m{row+1},'interpreter','none'); 
        ylabel(m{column+1},'interpreter','none');
        title('MUA');
        if row==1 && column==1
            legend({'','V1','V4'},'location','NorthWest');
        end
    end
end
if SaveFigs
    saveas(f6,fullfile(figfld, 'EPHYS_MUA_ModelComparison.png'));
end
if CloseFigs; close(f6); end

%% Stats model comparison R2 MUA
for RTH = [0 25]
    % stats V1
    MR2=R2(v1,2:5);
    sel=logical(sum(MR2>RTH,2));
    [p,tbl,stats] = kruskalwallis(MR2(sel,:),{'css','dog','p-lin','u-lin'});
    [c,m,h,gnames] = multcompare(stats);
    for i=1:size(c,1)
        fprintf(['V1 RTH-' num2str(RTH) ': ' gnames{c(i,1)} ' vs ' gnames{c(i,2)} ...
            ', p = ' num2str(c(i,6))  '\n'])
    end
    % stats V4
    MR2=R2(v4,2:5);
    sel=logical(sum(MR2>RTH,2));
    [p,tbl,stats] = kruskalwallis(MR2(sel,:),{'css','dog','p-lin','u-lin'});
    [c,m,h,gnames] = multcompare(stats);
    for i=1:size(c,1)
        fprintf(['V4 RTH-' num2str(RTH) ': ' gnames{c(i,1)} ' vs ' gnames{c(i,2)} ...
            ', p = ' num2str(c(i,6))  '\n'])
    end
end

%% LFP model comparison ===================================================
f7=figure; 
set(f7,'Position',[100 100 1600 1200]);
msz=15;

m=unique(tMUA.Model);
v1=tMUA.Area(strcmp(tMUA.Model,m{1}))==1;
v4=tMUA.Area(strcmp(tMUA.Model,m{1}))==4;

% scatter dots
sig=unique(tLFP.SigType);
lfp_order = [3 1 2 5 4];
spn=1; fbn=1;

for fb=lfp_order
    lfpmods{fbn,1}=[];
    lfpmods{fbn,2}=[];

    m=unique(tLFP.Model);
    % reorder ---------
    m=m([3 4 1 2]);
    modname={'P-LIN','U-LIN','CSS','DoG'};
    % -----------------
    R2=[];
    for i=1:length(m)
        R2 = [R2 tLFP.R2(...
            strcmp(tLFP.Model,m{i}) & ...
            strcmp(tLFP.SigType,sig{fb}))];
    end
        
    for m1=1:4
        lfpmods{fbn,1}=[lfpmods{fbn,1} R2(v1,m1)];
        lfpmods{fbn,2}=[lfpmods{fbn,2} R2(v4,m1)];
        for m2=m1+1:4
            subplot(length(sig),6,spn); hold on;
            plot([-5 100],[-5 100],'k');

%             M1R2=R2(v1,m1); 
%             M2R2=R2(v1,m2);           
%             M1R2=M1R2(M1R2>0 & M2R2>0);
%             M2R2=M2R2(M1R2>0 & M2R2>0);
            
            scatter(R2(v1,m1), R2(v1,m2),msz,'Marker','o',...
                'MarkerEdgeColor',[.3 .3 .3],'MarkerEdgeAlpha',0,...
                'MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',0.25);
%             scatter(R2(v4,m1), R2(v4,m2),msz,'Marker','o',...
%                 'MarkerEdgeColor',[.3 .8 .3],'MarkerEdgeAlpha',0,...
%                 'MarkerFaceColor',[.3 .8 .3],'MarkerFaceAlpha',0.25);

%             M1R2=R2(v4,m1); M1R2=M1R2(M1R2>0 & M2R2>0);
%             M2R2=R2(v4,m1); M1R2=M1R2(M1R2>0 & M2R2>0);          
                        
%             scatter(R2(v4,m1), R2(v4,m2),msz,'Marker','o',...
%                 'MarkerEdgeColor',[.3 .3 .3],'MarkerEdgeAlpha',0,...
%                 'MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',0.25);
            set(gca, 'Box','off', 'xlim', [-5 100], 'ylim',[-5 100],...
                'xticklabels',{},'yticklabels',{},'TickDir','out');
            xlabel(modname{m1},'interpreter','none');ylabel(modname{m2},'interpreter','none')
            title(sig{fb})
            if spn==1
                %legend({'','V1','V4'},'location','SouthEast');
            end
            spn=spn+1;
        end
    end
    fbn=fbn+1;
end
if SaveFigs
    saveas(f7,fullfile(figfld, 'EPHYS_LFP_ModelComparison.png'));
end
if CloseFigs; close(f7); end

%% stats ======
FreqNames = {'Theta','Alpha','Beta','Gamma-low','Gamma-high'};

for fb=1:5
    % stats V1
    fprintf(['\n= V1 LFP-' FreqNames{fb} ' ================\n'])
    [p,tbl,stats] = kruskalwallis(lfpmods{fb,1},modname);
    fprintf(['H =' num2str(tbl{2,5}(1)) ', df = ' num2str(tbl{2,3}(1)) ', p = ' num2str(tbl{2,6}(1)) '\n'])
    [c,m,h,gnames] = multcompare(stats);
    for i=1:size(c,1)
        fprintf([gnames{c(i,1)} ' vs ' gnames{c(i,2)} ...
            ', p = ' num2str(c(i,6))  '\n'])
    end
    % stats V4
    fprintf(['\n= V4 LFP-' FreqNames{fb} ' ================\n'])
    [p,tbl,stats] = kruskalwallis(lfpmods{fb,2},modname);
    fprintf(['H =' num2str(tbl{2,5}(1)) ', df = ' num2str(tbl{2,3}(1)) ', p = ' num2str(tbl{2,6}(1)) '\n'])
    [c,m,h,gnames] = multcompare(stats);
    for i=1:size(c,1)
        fprintf([gnames{c(i,1)} ' vs ' gnames{c(i,2)} ...
            ', p = ' num2str(c(i,6))  '\n'])
    end
end

%% R2 for different ephys signals =========================================
r2th=0;

ephys_MOD={'linear_ephys_cv1','linear_ephys_cv1_neggain',...
    'css_ephys_cv1','dog_ephys_cv1'};
ephys_MMS = MMS(:,1);

for m=1:2%1:length(ephys_MOD)
    RR=[];

    fprintf(['\n============ ' ephys_MOD{m} ' ===========\n'])
    
    model=ephys_MOD{m};
    s = strcmp(tMUA.Model,model);
    RR=[RR tMUA.R2(s)];
    
    s = strcmp(tLFP.Model,model);
    sig=unique(tLFP.SigType);
    lfp_order = [3 1 2 5 4];
    for i=lfp_order
        b=strcmp(tLFP.SigType,sig{i});
        RR=[RR tLFP.R2(s & b)];
    end
    LAB=['MUA';sig(lfp_order)];
    
    f8=figure; 
    set(f8,'Position',[100 100 1900 1350]);
    sgtitle(['R2 per Model: ' model],'interpreter','none');
    
    c=0;d=0;
    for ref=1:6
        c=c+1;
        for fb=1:6
            d=d+1;
            s=(RR(:,ref)>r2th & RR(:,fb)>r2th);
            subplot(6,6,d); hold on; 
            %scatter(RR(s,ref),RR(s,fb),120,[0.3 0.3 0.3],'Marker','.');
            binscatter(RR(s,ref),RR(s,fb),25,...
                'XLimits', [0 100],...
                'YLimits', [0 100],...
                'ShowEmptyBins', 'off')
            colorbar; 
            set(gca,'ColorScale','log');
            colormap(inferno)
            caxis([1 256]) 
            plot([0 100],[0 100],'Color',[.7 .7 .7],'LineWidth',2);
            xlabel(LAB{ref});ylabel(LAB{fb});
            %title(['R2 ' model],'interpreter','none');
            set(gca,'xlim',[0 100],'ylim',[0 100]);
            set(gca,'TickDir','out','xtick',[0 50 100],'ytick',[0 50 100]);
        end
    end
    if SaveFigs
        saveas(f8,fullfile(figfld, ['EPHYS_MUA_R2_' ephys_MOD{m} '.png']));
    end
    if CloseFigs; close(f8); end
    
%     % Distance from diagonal ==============================================
%     f9=figure; set(f9,'Position',[100 100 1300 1200]);
%     LAB=['MUA';sig(lfp_order)];
%     sgtitle(['Differences Model: ' model],'interpreter','none');
% 
%     c=0;d=0;
%     for ref=1:6
%         c=c+1;
%         for fb=1:6
%             d=d+1;
%             s=(RR(:,ref)>r2th & RR(:,fb)>r2th);
%             
%             subplot(6,6,d); hold on;
%             dRR = RR(s,fb)-RR(s,ref);
%             h = histogram(dRR,-100:1:100,'FaceColor','k','FaceAlpha',1);
%             YY = get(gca,'ylim');
%             plot([0 0],YY,'Color',[0.5 0.5 0.5],'LineWidth',2)
%             plot([mean(dRR) mean(dRR)],YY,'r','LineWidth',2)
%             plot([median(dRR) median(dRR)],YY,'b','LineWidth',2)
%             set(gca,'xlim',[-100 100])
%             xlabel(['dRR ' model],'interpreter','none');
%             ylabel('cnt','interpreter','none');
%             title(['R2 ' LAB{(fb)} '-' LAB{(ref)} ],...
%                 'interpreter','none')
%             
%             if ref ==1 && fb ==1
%                 legend({'HIST','0','MEAN','MEDIAN'});
%             end
%             
%         end
%     end
%     if SaveFigs
%         saveas(f9,fullfile(figfld, ['EPHYS_MUA_R2diff_' ephys_MOD{m} '.png']));  
%     end
%     if CloseFigs; close(f9); end


    % Stats
    [p,tbl,stats] = kruskalwallis(RR,LAB');
    fprintf(['H =' num2str(tbl{2,5}(1)) ', df = ' num2str(tbl{2,3}(1)) ', p = ' num2str(tbl{2,6}(1)) '\n'])
    [c,m,h,gnames] = multcompare(stats);
    for i=1:size(c,1)
        fprintf([gnames{c(i,1)} ' (' num2str(m(c(i,1),1)) ' +/- ' num2str(m(c(i,1),2)) ') vs ' ...
            gnames{c(i,2)} ' (' num2str(m(c(i,2),1)) ' +/- ' num2str(m(c(i,2),2)) '), p = ' num2str(c(i,6))  '\n'])
    end    
    close all
end

%% pRF size for different ephys signals ===================================
% SZ is [ MUA_R2(1) MUA_RFS(2) 
%         THETA_R2(3) THETA_RFS(4) 
%         ALPHA_R2(5) ALPHA_RFS(6) 
%         BETA_R2(7) BETA_RFS8) 
%         LGAM_R2(9) LG_RFS(10) 
%         HGAM_R2(11) HGAM_RFS(12)]
r2th=25;

ephys_MOD={'linear_ephys_cv1','linear_ephys_cv1_neggain',...
    'css_ephys_cv1','dog_ephys_cv1'};
ephys_MMS = MMS(:,1);

for m=1:length(ephys_MOD)
    model=ephys_MOD{m};
    LAB=['MUA';sig(lfp_order)];
    
    f10=figure; set(f10,'Position',[100 100 1300 1200]);
    sgtitle(['SZ per Model: ' model],'interpreter','none');
    
    c=0;d=0;
    for ref=1:2:12
        c=c+1;
        SS(c).nSZ =[];
        for fb=1:2:12
            d=d+1;
            s=(SZ{m}(:,ref)>r2th & SZ{m}(:,fb)>r2th);
            
            subplot(6,6,d); hold on; plot([0 100],[0 100],'k');
            
            scatter(SZ{m}(s,ref+1),SZ{m}(s,fb+1),120,[0.3 0.3 0.3],'Marker','.')
            xlabel(LAB{(ref+1)/2});ylabel(LAB{(fb+1)/2});
            title(['pRF size ' model(1:3)],'interpreter','none')
            
            set(gca,'xlim',[0 10],'ylim',[0 10]);
            
            SS(c).nSZ = [SS(c).nSZ ; ...
                median(diff(SZ{m}(s,[ref+1 fb+1]),1,2)) ...
                std(diff(SZ{m}(s,[ref+1 fb+1]),1,2))./sqrt(sum(s)) ...
                median(  SZ{m}(s,fb+1)./SZ{m}(s,ref+1) ) ...
                std(  SZ{m}(s,fb+1)./SZ{m}(s,ref+1) )./sqrt(sum(s))];
        end
    end
    if SaveFigs
        saveas(f10,fullfile(figfld, ['EPHYS_SZ_' ephys_MOD{m} '.png']));  
    end
    if CloseFigs; close(f10); end
    
    % Distance from diagonal ==============================================
    
    f11=figure; set(f11,'Position',[100 100 1300 1200]);
    sgtitle(['SZ DIFF per Model: ' model],'interpreter','none');
    
    c=0;d=0;
    for ref=1:2:12
        c=c+1;
        SS(c).nSZ =[];
        for fb=1:2:12
            d=d+1;
            s=(SZ{m}(:,ref)>r2th & SZ{m}(:,fb)>r2th);
            
            subplot(6,6,d); hold on;
            dSz = SZ{m}(s,fb+1)-SZ{m}(s,ref+1);
            h = histogram(dSz,-5:0.1:5,'FaceColor','k','FaceAlpha',1);
            YY = get(gca,'ylim');
            plot([0 0],YY,'Color',[0.5 0.5 0.5],'LineWidth',2)
            plot([mean(dSz) mean(dSz)],YY,'r','LineWidth',2)
            plot([median(dSz) median(dSz)],YY,'b','LineWidth',2)
            set(gca,'xlim',[-6 6])
            xlabel(['Sz Diff - ' model(1:3)],'interpreter','none');
            ylabel('cnt','interpreter','none');
            title(['Size ' LAB{(fb+1)/2} '-' LAB{(ref+1)/2} ],...
                'interpreter','none')
            
            if ref ==1 && fb ==1
                legend({'HIST','0','MEAN','MEDIAN'});
            end
            
        end
    end
    if SaveFigs
        saveas(f11,fullfile(figfld, ['EPHYS_SZ-diff_' ephys_MOD{m} '.png']));
    end
    if CloseFigs; close(f11); end
end





%% ========================================================================
%  FMRI-EPHYS PRF ANALYSIS ------------------------------------------------
%  ========================================================================

%% Correlate MRI-ephys RESAMPLE ON MAP-GRID ===============================
% based on the full map X,Y,size
% This analysis takes a while!! Do not overuse...
rng(1); % seed the random number generator

Rth_mri = 5; % R2 threshold MRI
Rth_ephys = 30; % R2 threshold ephys
mxS = 10; % maximum size

MODS = {...
    'linhrf_cv1_mhrf','linear_ephys_cv1';...
    'linhrf_cv1_mhrf_neggain','linear_ephys_cv1_neggain';...
    'csshrf_cv1_mhrf','css_ephys_cv1';...
    'doghrf_cv1_mhrf','dog_ephys_cv1';...
    };
MRI_MODEL = MODS(:,1);
EPHYS_MODEL = MODS(:,2);
MMS={'linear','linear_ng','css','dog'};

nbtstr = 200;
np = 100;
grid_vf = [ 0 5 -5 0 ; 0 8 -8 0]; % [xmin xmax ymin ymax] [v1; v4] dva
grid_spacing = 0.25;% dva
pth = 01;
poscorr_only = true;

warning off;

cmROI = {'V1','V4'};
fprintf('=======================\n');
for m = 1:size(MODS,1)
    fprintf(['\nBootstrap Correlation for Model: ' MODS{m} '\n']);
    
    s_R2 = T(modidx.(MRI_MODEL{m})).mod.R2 > Rth_mri & ...
         T(modidx.(MRI_MODEL{m})).mod.rfs < mxS;
    
    % collect mri prfs
    for r = 1:size(cmROI,2)
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{m})).mod.ROI,...
                ck_GetROIidx(cmROI(r),rois) );
        if strcmp(cmROI{r},'V1') % V1
            mri1(m).X = T(modidx.(MRI_MODEL{m})).mod.X(SSS);
            mri1(m).Y = T(modidx.(MRI_MODEL{m})).mod.Y(SSS);
            if strcmp(MRI_MODEL{m}(1:3),'dog')
                mri1(m).S = T(modidx.(MRI_MODEL{m})).mod.rfs(SSS);
            else
                mri1(m).S = T(modidx.(MRI_MODEL{m})).mod.rfs(SSS);
            end
        elseif strcmp(cmROI{r},'V4') % V4
            mri4(m).X = T(modidx.(MRI_MODEL{m})).mod.X(SSS);
            mri4(m).Y = T(modidx.(MRI_MODEL{m})).mod.Y(SSS);
            mri4(m).S = T(modidx.(MRI_MODEL{m})).mod.rfs(SSS);
        end
    end
               
    % collect ephys prfs
    % MUA V1
    s = strcmp(tMUA.Model,EPHYS_MODEL{m}) & ...
        tMUA.Area == 1 & tMUA.R2 > Rth_ephys & tMUA.rfs < mxS;
    mua1(m).X = tMUA.X(s); 
    mua1(m).Y = tMUA.Y(s); 
    mua1(m).S = tMUA.rfs(s);
    % MUA V4
    s = strcmp(tMUA.Model,EPHYS_MODEL{m}) & ...
        tMUA.Area == 4 & tMUA.R2 > Rth_ephys & tMUA.rfs < mxS;
    mua4(m).X = tMUA.X(s); 
    mua4(m).Y = tMUA.Y(s); 
    mua4(m).S = tMUA.rfs(s);  
    % LFP
    freqband=unique(tLFP.SigType);
    for fb = 1: length(freqband)
        % V1
        s = strcmp(tLFP.Model,EPHYS_MODEL{m}) & ...
            tLFP.Area == 1 & tLFP.R2 > Rth_ephys & tLFP.rfs < mxS & ...
            strcmp(tLFP.SigType, freqband{fb});
        lfp1(fb,m).freqband =  freqband{fb};
        lfp1(fb,m).X =  tLFP.X(s);
        lfp1(fb,m).Y =  tLFP.Y(s);
        lfp1(fb,m).S =  tLFP.rfs(s);
        
        % V4
        s = strcmp(tLFP.Model,EPHYS_MODEL{m}) & ...
            tLFP.Area == 4 & tLFP.R2 > Rth_ephys & tLFP.rfs < mxS & ...
            strcmp(tLFP.SigType, freqband{fb});
        lfp4(fb,m).freqband =  freqband{fb};
        lfp4(fb,m).X =  tLFP.X(s);
        lfp4(fb,m).Y =  tLFP.Y(s);
        lfp4(fb,m).S =  tLFP.rfs(s);
    end
   
    % Calculate  grids and do bootstrap correlation
    % =============================================
    % create XY-grid
    [x1q,y1q] = meshgrid(...
        grid_vf(1,1):grid_spacing:grid_vf(1,2),...
        grid_vf(1,3):grid_spacing:grid_vf(1,4));
    [x4q,y4q] = meshgrid(...
        grid_vf(2,1):grid_spacing:grid_vf(2,2),...
        grid_vf(2,3):grid_spacing:grid_vf(2,4));
    
    mri1(m).S_grid = griddata(mri1(m).X,mri1(m).Y,mri1(m).S,x1q,y1q,'linear');
    mri4(m).S_grid = griddata(mri4(m).X,mri4(m).Y,mri4(m).S,x4q,y4q,'linear');
    mua1(m).S_grid = griddata(mua1(m).X,mua1(m).Y,mua1(m).S,x1q,y1q,'linear');
    mua4(m).S_grid = griddata(mua4(m).X,mua4(m).Y,mua4(m).S,x4q,y4q,'linear');
    for fb = 1: length(freqband)
        lfp1(fb,m).S_grid = griddata(...
            lfp1(fb,m).X,lfp1(fb,m).Y,lfp1(fb,m).S,x1q,y1q,'linear');
        lfp4(fb,m).S_grid = griddata(...
            lfp4(fb,m).X,lfp4(fb,m).Y,lfp4(fb,m).S,x4q,y4q,'linear');
    end
    
    if false % true
        figure;
        subplot(2,2,1); hold on;
        contourf(x1q,y1q,mri1(m).S_grid,'LevelStep',0.1,'LineStyle','none');
        plot(x1q,y1q,'k.')
        
        subplot(2,2,2);hold on;
        contourf(x1q,y1q,mua1(m).S_grid,'LevelStep',0.1,'LineStyle','none');
        plot(x1q,y1q,'k.')
        
        subplot(2,2,3);hold on;
        contourf(x4q,y4q,mri4(m).S_grid,'LevelStep',0.1,'LineStyle','none');
        plot(x4q,y4q,'k.')
        
        subplot(2,2,4);hold on;
        contourf(x4q,y4q,mua4(m).S_grid,'LevelStep',0.1,'LineStyle','none');
        plot(x4q,y4q,'k.')
    end
        
    % bootstrap the correlation analysis
    fprintf('Performing a bootstrapped correlation analysis\n')
    v1 = find(~isnan(mri1(m).S_grid));
    v4 = find(~isnan(mri4(m).S_grid));
    cc1=[]; cc4=[];
    pp1=[]; pp4=[];
    
    fprintf(['nBtstr: ' num2str(nbtstr) ', nSamples: ' num2str(np) '\n'])
    plotscatter=false;
    if plotscatter; figure; end % plotting all here, selecting later
    for i=1:nbtstr
        % --- V1 ---
        c1=[];p1=[];
        V=v1(randperm(length(v1)));
        selS = [mri1(m).S_grid(:) mua1(m).S_grid(:)];
        selS = selS(V(1:np),:);
        selS = selS(~isnan(selS(:,2)),:);
        [r,p]=corrcoef(selS(:,1), selS(:,2));
        if plotscatter
            subplot(2,6,1); hold on;
            scatter(mri1(m).S_grid(V(1:np)), mua1(m).S_grid(V(1:np)),'o');
            title('V1 Map corr.');xlabel('MRI');ylabel('MUA');
            set(gca,'xlim',[0 6],'ylim',[0 6]);
        end
        c1=[c1 r(2)]; p1=[p1 p(2)];
        
        for fb=1:length(freqband)
            try
                selS = [mri1(m).S_grid(:) lfp1(fb,m).S_grid(:)];
                selS = selS(V(1:np),:);
                selS = selS(~isnan(selS(:,2)),:);
                [r,p]=corrcoef(selS(:,1), selS(:,2));
                if plotscatter
                    subplot(2,6,1+fb);hold on;
                    scatter(mri1(m).S_grid(V(1:np)), lfp1(fb,m).S_grid(V(1:np)),'o');
                    set(gca,'xlim',[0 6],'ylim',[0 6]);
                end
            catch
                r=NaN(2,2);
            end
            if plotscatter
                subplot(2,6,1+fb);
                title('V1 Map corr.');xlabel('MRI');ylabel(lfp1(fb,m).freqband);
            end
            c1=[c1 r(2)]; p1=[p1 p(2)];
        end
        cc1=[cc1; c1]; pp1=[pp1; p1];
        stats(m).cc1 = cc1; stats(m).pp1 = pp1;
        
        % --- V4 ----
        c4=[];p4=[];
        V=v4(randperm(length(v4)));
        
        selS = [mri4(m).S_grid(:) mua4(m).S_grid(:)];
        selS = selS(V(1:np),:);
        selS = selS(~isnan(selS(:,2)),:);
        [r,p]=corrcoef(selS(:,1), selS(:,2));
        if plotscatter
            subplot(2,6,7);hold on;
            scatter(mri4(m).S_grid(V(1:np)), mua4(m).S_grid(V(1:np)),'o');
            title('V4 Map corr.');xlabel('MRI');ylabel('MUA');
            set(gca,'xlim',[0 6],'ylim',[0 6]);
        end
        c4=[c4 r(2)]; p4=[p4 p(2)];
        
        for fb=1:length(freqband)
            try
                selS = [mri4(m).S_grid(:) lfp4(fb,m).S_grid(:)];
                selS = selS(V(1:np),:);
                selS = selS(~isnan(selS(:,2)),:);
                [r,p]=corrcoef(selS(:,1), selS(:,2));
                if plotscatter
                    subplot(2,6,7+fb);hold on;
                    scatter(mri4(m).S_grid(V(1:np)), lfp4(fb,m).S_grid(V(1:np)),'o');
                    set(gca,'xlim',[0 6],'ylim',[0 6]);
                end
            catch
                r=NaN(2,2);
            end
            if plotscatter
                subplot(2,6,7+fb);
                title('V4 Map corr.');xlabel('MRI');ylabel(lfp1(fb,m).freqband);
            end
            if length(p)>1
                c4=[c4 r(2)]; p4=[p4 p(2)];
            else
                 c4=[c4 nan]; p4=[p4 nan];
            end
        end
        cc4=[cc4; c4]; pp4=[pp4; p4];
        stats(m).cc4 = cc4; stats(m).pp4 = pp4;
    end
    
    stats(m).c1filt=stats(m).cc1; stats(m).c4filt=stats(m).cc4;
    if poscorr_only
        stats(m).c1filt(stats(m).pp1>pth | stats(m).cc1<0)=NaN;
        stats(m).c4filt(stats(m).pp4>pth | stats(m).cc4<0)=NaN;
    else
        stats(m).c1filt(stats(m).pp1>pth)=NaN;
        stats(m).c4filt(stats(m).pp4>pth)=NaN;
    end
    stats(m).cc1_stat = [nanmean(stats(m).c1filt,1); nanstd(stats(m).c1filt,0,1)];
    stats(m).cc4_stat = [nanmean(stats(m).c4filt,1); nanstd(stats(m).c4filt,0,1)];
    
    % plot average and stdev
    f14=figure; hold on;
    hBar = bar(1:6, [stats(m).cc1_stat(1,:);stats(m).cc4_stat(1,:)]);pause(0.1)
    xBar=cell2mat(get(hBar,'XData')).' + [hBar.XOffset];
    errorbar(xBar, [stats(m).cc1_stat(1,:);stats(m).cc4_stat(1,:)]',...
        [stats(m).cc1_stat(2,:);stats(m).cc4_stat(2,:)]','k','LineStyle','none')
    GroupLabels = {'MUA',...
        lfp1(1,m).freqband,...
        lfp1(2,m).freqband,...
        lfp1(3,m).freqband,...
        lfp1(4,m).freqband,...
        lfp1(5,m).freqband};
    set(gca,'xtick',1:6,'xticklabels',GroupLabels);
    legend({'V1','V4'},'Location','NorthWest');
    title(...
        {['pRF model: ' MMS{m}],...
        ['nPoints: ' num2str(np) ', nBtstr: ' num2str(nbtstr) ...
        ', p < ' num2str(pth) ],[' R2th-mri: '  num2str(Rth_mri) ...
        ', R2th-ephys: ' num2str(Rth_ephys)]},'Interpreter','none')
    
    if SaveFigs 
        saveas(f14,fullfile(figfld, ['XMOD_GRID_' ephys_MOD{m} '.png']));  
    end
    
    % Stats ---
    if false
        % V1
        [stats(m).p1,stats(m).t1,stats(m).s1] = ...
            anova1(stats(m).c1filt,GroupLabels);
        [stats(m).comp1, stats(m).means1, stats(m).h1, stats(m).names1] = ...
            multcompare(stats(m).s1);
        % V4
        [stats(m).p4,stats(m).t4,stats(m).s4] = ...
            anova1(stats(m).c4filt,GroupLabels);
        [stats(m).comp4, stats(m).means4, stats(m).h4, stats(m).names4] = ...
            multcompare(stats(m).s4);
    end
    if CloseFigs; close(f14); end
end
warning on;

%% Correlate MRI-ephys REGRESSION BASED ===================================
% based on the ecc vs size correlation
% PER SIGNAL SOURCE
% - calculate linear regression on random (filtered) subset of data
% - repeat
% - filter results
% ACROSS SIGNAL SOURCES
% - correlate the bootstrapped fit-result parameters
% - optional >> do this in a bootstrapped way

% This analysis takes a while!! Do not overuse...
rng(1); % seed the random number generator

Rth_mri = 10; % R2 threshold MRI
Rth_ephys = 25; % R2 threshold ephys
mxS = 1000;%25; % maximum size
MaxECC = 25; % max ecc to use for fitting
selVFC = true;

MODS = {...
    'linhrf_cv1_mhrf','linear_ephys_cv1';...
    'linhrf_cv1_mhrf_neggain','linear_ephys_cv1_neggain';...
    'csshrf_cv1_mhrf','css_ephys_cv1';...
    'doghrf_cv1_mhrf','dog_ephys_cv1';...
    };
MRI_MODEL = MODS(:,1); EPHYS_MODEL = MODS(:,2);
MMS={'linear','linear_ng','css','dog'};

nbtstr = 200;
np = 100;
pth = 01;
poscorr_only = false;

warning off;
cmROI = {'V1','V4'};
fprintf('=======================\n');
for m = 3%[1 2 3] % 1:size(MODS,1)
    fprintf(['\nCrossmodal Correlation for Model: ' MODS{m} '\n']);
    
    s_R2 = T(modidx.(MRI_MODEL{m})).mod.R2 > Rth_mri & ...
        T(modidx.(MRI_MODEL{m})).mod.rfs < mxS;
    
    % collect mri prfs
    for r = 1:size(cmROI,2)
        SSS = s_R2 & ismember( T(modidx.(MRI_MODELS{m})).mod.ROI,...
            ck_GetROIidx(cmROI(r),rois) );
        
        if selVFC
            VFSEL = T(modidx.(MRI_MODEL{m})).mod.X >= 0 & ...
                T(modidx.(MRI_MODEL{m})).mod.Y <=0;
            SSS = SSS & VFSEL;
        end
        
        if strcmp(cmROI{r},'V1') % V1
            mri1(m).ECC = T(modidx.(MRI_MODEL{m})).mod.ecc(SSS);
            mri1(m).S = T(modidx.(MRI_MODEL{m})).mod.rfs(SSS);
            mri1(m).G = T(modidx.(MRI_MODEL{m})).mod.gain(SSS);
        elseif strcmp(cmROI{r},'V4') % V4
            mri4(m).ECC = T(modidx.(MRI_MODEL{m})).mod.ecc(SSS);
            mri4(m).S = T(modidx.(MRI_MODEL{m})).mod.rfs(SSS);
            mri4(m).G = T(modidx.(MRI_MODEL{m})).mod.gain(SSS);
        end
    end
    
    % collect ephys prfs
    % MUA V1
    s = strcmp(tMUA.Model,EPHYS_MODEL{m}) & ...
        tMUA.Area == 1 & tMUA.R2 > Rth_ephys & tMUA.rfs < mxS;
    % exclude the outlier V1 arrays in both  monkeys (they might be in WM)
%     s2 = ( strcmp(tMUA.Monkey,'lick') & tMUA.Array == 11) | ...
%         ( strcmp(tMUA.Monkey,'aston') & tMUA.Array == 10);
    % Exclude 2 outlier V1 arrays in both monkeys
    s2 = ( strcmp(tMUA.Monkey,'lick') & (tMUA.Array == 11 | tMUA.Array == 13)) | ...
        ( strcmp(tMUA.Monkey,'aston') & (tMUA.Array == 10 | tMUA.Array == 12));
    s(s2) = false;
    
    mua1(m).ECC = tMUA.ecc(s);
    mua1(m).S = tMUA.rfs(s);
    mua1(m).G = tMUA.gain(s);
    % MUA V4
    s = strcmp(tMUA.Model,EPHYS_MODEL{m}) & ...
        tMUA.Area == 4 & tMUA.R2 > Rth_ephys & tMUA.rfs < mxS;
    mua4(m).ECC = tMUA.ecc(s);
    mua4(m).S = tMUA.rfs(s);
    mua4(m).G = tMUA.gain(s);
    
    % LFP
    freqband=unique(tLFP.SigType);
    for fb = 1: length(freqband)
        % V1
        s = strcmp(tLFP.Model,EPHYS_MODEL{m}) & ...
            tLFP.Area == 1 & tLFP.R2 > Rth_ephys & tLFP.rfs < mxS & ...
            strcmp(tLFP.SigType, freqband{fb});
%         % exclude the outlier V1 arrays in both  monkeys (they might be in WM)
%         s2 = ( strcmp(tLFP.Monkey,'lick') & tLFP.Array == 11) | ...
%             ( strcmp(tLFP.Monkey,'aston') & tLFP.Array == 10);
        % Exclude 2 outlier V1 arrays in both monkeys
        s2 = ( strcmp(tMUA.Monkey,'lick') & (tMUA.Array == 11 | tMUA.Array == 13)) | ...
            ( strcmp(tMUA.Monkey,'aston') & (tMUA.Array == 10 | tMUA.Array == 12));
        s(s2) = false;
        
        lfp1(fb,m).freqband =  freqband{fb};
        lfp1(fb,m).ECC =  tLFP.ecc(s);
        lfp1(fb,m).S =  tLFP.rfs(s);
        lfp1(fb,m).G =  tLFP.gain(s);
        
        % V4
        s = strcmp(tLFP.Model,EPHYS_MODEL{m}) & ...
            tLFP.Area == 4 & tLFP.R2 > Rth_ephys & tLFP.rfs < mxS & ...
            strcmp(tLFP.SigType, freqband{fb});
        lfp4(fb,m).freqband =  freqband{fb};
        lfp4(fb,m).ECC =  tLFP.ecc(s);
        lfp4(fb,m).S =  tLFP.rfs(s);
        lfp4(fb,m).G =  tLFP.gain(s);
    end
    
    % Calculate  & bootstrap linear regressions
    % =============================================
    fprintf('Performing regressions on ALL data\n')
    MRI_stats1 = regstats(mri1(m).S,mri1(m).ECC);
    MRI_stats4 = regstats(mri4(m).S,mri4(m).ECC);
    MUA_stats1 = regstats(mua1(m).S,mua1(m).ECC);
    MUA_stats4 = regstats(mua4(m).S,mua4(m).ECC);
    
    for fb=1:length(freqband)
        try
            LFP_stats1{fb} = regstats(lfp1(fb,m).S,lfp1(fb,m).ECC);
        catch
            LFP_stats1{fb} = [];
        end
        try
            LFP_stats4{fb} = regstats(lfp4(fb,m).S,lfp4(fb,m).ECC);
        catch
            LFP_stats4{fb} = [];
        end
    end
    
    % Do statistics on model comparisons
    % use fitlm
    
    ck_xmod_stats;
    XMOD(m).model = MMS{m}; 
    XMOD(m).stats = stats;
    
    if m==2
        ck_xmod_stats_lowfreqlfp;
    end
    
    % =====================================================================
    DoBootstrapApproach = false;
    if DoBootstrapApproach
        fprintf('Performing a bootstrapped regression analysis\n')
        cc1=[]; cc4=[];
        pp1=[]; pp4=[];
        
        fprintf(['nBtstr: ' num2str(nbtstr) ', nSamples: ' num2str(np) '\n'])
        for i=1:nbtstr
            c1=[];p1=[];
            % --- V1 ---
            % mri
            V=randperm(length(mri1(m).ECC));
            selS = V(1:np);
            stats = regstats(mri1(m).S(selS),mri1(m).ECC(selS));
            c1=[c1 stats.beta'];
            p1=[p1 stats.tstat.pval'];
            
            % mua
            V=randperm(length(mua1(m).ECC));
            selS = V(1:np);
            stats = regstats(mua1(m).S(selS),mua1(m).ECC(selS));
            c1=[c1 stats.beta'];
            p1=[p1 stats.tstat.pval'];
            
            for fb=1:length(freqband)
                try
                    V=randperm(length(lfp1(fb,m).ECC));
                    selS = V(1:np);
                    stats = regstats(lfp1(fb,m).S(selS),lfp1(fb,m).ECC(selS));
                    c1=[c1 stats.beta'];
                    p1=[p1 stats.tstat.pval'];
                catch
                    c1=[c1 nan nan];
                    p1=[p1 nan nan];
                    %fprintf(['fb' num2str(fb) ': error with regstats\n'])
                end
            end
            cc1=[cc1; c1]; pp1=[pp1; p1];
            
            % --- V4 ----
            c4=[];p4=[];
            % mri
            V=randperm(length(mri4(m).ECC));
            selS = V(1:np);
            stats = regstats(mri4(m).S(selS),mri4(m).ECC(selS));
            c4=[c4 stats.beta'];
            p4=[p4 stats.tstat.pval'];
            
            % mua
            V=randperm(length(mua4(m).ECC));
            selS = V(1:np);
            stats = regstats(mua4(m).S(selS),mua4(m).ECC(selS));
            c4=[c4 stats.beta'];
            p4=[p4 stats.tstat.pval'];
            
            for fb=1:length(freqband)
                try
                    V=randperm(length(lfp4(fb,m).ECC));
                    selS = V(1:np);
                    stats = regstats(lfp4(fb,m).S(selS),lfp4(fb,m).ECC(selS));
                    c4=[c4 stats.beta'];
                    p4=[p4 stats.tstat.pval'];
                catch
                    c4=[c4 nan nan];
                    p4=[p4 nan nan];
                    %fprintf(['fb' num2str(fb) ': error with regstats\n'])
                end
            end
            cc4=[cc4; c4]; pp4=[pp4; p4];
        end
        
        % show slope distributions for all signals
        pthr=0.05;
        
        figure;
        subplot(2,4,1)
        histogram(cc1(pp1(:,2)<pthr,2),-0.2:0.02:0.5)
        title('MRI V1'); xlabel('Ecc-Sz slope')
        subplot(2,4,2)
        histogram(cc1(pp1(:,4)<pthr,4),-0.2:0.02:0.5)
        title('MUA V1');xlabel('Ecc-Sz slope')
        subplot(2,4,3)
        histogram(cc1(pp1(:,6)<pthr,6),-0.2:0.02:0.5)
        title('ALPHA V1');xlabel('Ecc-Sz slope')
        subplot(2,4,4)
        histogram(cc1(pp1(:,8)<pthr,8),-0.2:0.02:0.5)
        title('BETA V1');xlabel('Ecc-Sz slope')
        subplot(2,4,5)
        histogram(cc1(pp1(:,10)<pthr,10),-0.2:0.02:0.5)
        title('THETA V1');xlabel('Ecc-Sz slope')
        subplot(2,4,6)
        histogram(cc1(pp1(:,12)<pthr,12),-0.2:0.02:0.5)
        title('H-GAMMA V1');xlabel('Ecc-Sz slope')
        subplot(2,4,7)
        histogram(cc1(pp1(:,14)<pthr,14),-0.2:0.02:0.5)
        title('L-GAMMA V1');xlabel('Ecc-Sz slope')
        subplot(2,4,8)
        ae_cc1 = [cc1(:,1:2); cc1(:,3:4); cc1(:,5:6); cc1(:,7:8);...
            cc1(:,9:10); cc1(:,11:12); cc1(:,13:4)];
        ae_pp1 = [pp1(:,1:2); pp1(:,3:4); pp1(:,5:6); pp1(:,7:8);...
            pp1(:,9:10); pp1(:,11:12); pp1(:,13:4)];
        histogram(ae_cc1(ae_pp1(:,2)<pthr,2),-0.2:0.02:0.5)
        title('ALL EPHYS V1');xlabel('Ecc-Sz slope')
        sgtitle(['MODELS: ' MODS{m,1} ' and ' MODS{m,2}],'interpreter','none')
        
        
        figure;
        subplot(2,4,1)
        histogram(cc4(pp4(:,2)<pthr,2),-0.2:0.02:0.5)
        title('MRI V4');xlabel('Ecc-Sz slope')
        subplot(2,4,2)
        histogram(cc4(pp4(:,4)<pthr,4),-0.2:0.02:0.5)
        title('MUA V4');xlabel('Ecc-Sz slope')
        subplot(2,4,3)
        histogram(cc4(pp4(:,6)<pthr,6),-0.2:0.02:0.5)
        title('ALPHA V4');xlabel('Ecc-Sz slope')
        subplot(2,4,4)
        histogram(cc4(pp4(:,8)<pthr,8),-0.2:0.02:0.5)
        title('BETA V4');xlabel('Ecc-Sz slope')
        subplot(2,4,5)
        histogram(cc4(pp4(:,10)<pthr,10),-0.2:0.02:0.5)
        title('THETA V4');xlabel('Ecc-Sz slope')
        subplot(2,4,6)
        histogram(cc4(pp4(:,12)<pthr,12),-0.2:0.02:0.5)
        title('H-GAMMA V4');xlabel('Ecc-Sz slope')
        subplot(2,4,7)
        histogram(cc4(pp4(:,14)<pthr,14),-0.2:0.02:0.5)
        title('L-GAMMA V4');xlabel('Ecc-Sz slope')
        subplot(2,4,8)
        ae_cc4 = [cc4(:,1:2); cc4(:,3:4); cc4(:,5:6); cc4(:,7:8);...
            cc4(:,9:10); cc4(:,11:12); cc4(:,13:4)];
        ae_pp4 = [pp4(:,1:2); pp4(:,3:4); pp1(:,5:6); pp4(:,7:8);...
            pp4(:,9:10); pp4(:,11:12); pp4(:,13:4)];
        histogram(ae_cc4(ae_pp4(:,2)<pthr,2),-0.2:0.02:0.5)
        title('ALL EPHYS V4');xlabel('Ecc-Sz slope')
        sgtitle(['MODELS: ' MODS{m,1} ' and ' MODS{m,2}],'interpreter','none')
        
        
        % plot average and stdev
        f14=figure; hold on;
        hBar = bar(1:6, [stats(m).cc1_stat(1,:);stats(m).cc4_stat(1,:)]);pause(0.1)
        xBar=cell2mat(get(hBar,'XData')).' + [hBar.XOffset];
        errorbar(xBar, [stats(m).cc1_stat(1,:);stats(m).cc4_stat(1,:)]',...
            [stats(m).cc1_stat(2,:);stats(m).cc4_stat(2,:)]','k','LineStyle','none')
        GroupLabels = {'MUA',...
            lfp1(1,m).freqband,...
            lfp1(2,m).freqband,...
            lfp1(3,m).freqband,...
            lfp1(4,m).freqband,...
            lfp1(5,m).freqband};
        set(gca,'xtick',1:6,'xticklabels',GroupLabels);
        legend({'V1','V4'},'Location','NorthWest');
        title(...
            {['pRF model: ' MMS{m}],...
            ['nPoints: ' num2str(np) ', nBtstr: ' num2str(nbtstr) ...
            ', p < ' num2str(pth) ],[' R2th-mri: '  num2str(Rth_mri) ...
            ', R2th-ephys: ' num2str(Rth_ephys)]},'Interpreter','none')
        
        %saveas(f14,fullfile(figfld, ['MRI-EPHYS_' ephys_MOD{m} '.png']));
        
        % Stats ---
        if false
            % V1
            [stats(m).p1,stats(m).t1,stats(m).s1] = ...
                anova1(stats(m).c1filt,GroupLabels);
            [stats(m).comp1, stats(m).means1, stats(m).h1, stats(m).names1] = ...
                multcompare(stats(m).s1);
            % V4
            [stats(m).p4,stats(m).t4,stats(m).s4] = ...
                anova1(stats(m).c4filt,GroupLabels);
            [stats(m).comp4, stats(m).means4, stats(m).h4, stats(m).names4] = ...
                multcompare(stats(m).s4);
        end
    end
end
warning on;

%% What's specific about the good DoG fits MUA EDITION ====================
% - find these channels
% - plot their location
% - plot their size
% - plot their prf profile

R2th = 25; % minimum R2
R2enh = 5; % R2 improvement
   
DoG = tMUA(...
    strcmp(tMUA.Model,'dog_ephys_cv1'),:);
lin_n = tMUA(...
    strcmp(tMUA.Model,'linear_ephys_cv1_neggain'),:);
lin = tMUA(...
    strcmp(tMUA.Model,'linear_ephys_cv1'),:);

f_neg1 = figure;
set(f_neg1,'Position',[10 10 1200 500]);
chan_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh;
subplot(1,2,1);scatter(DoG.X(chan_sel),DoG.Y(chan_sel),...
    'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
set(gca,'xaxislocation','origin','yaxislocation','origin',...
    'xlim',[-8 8],'ylim',[-8 8]);
title('Locations of pRF with good DoG fits & bad LIN fits')
xlabel('X deg');ylabel('Y deg');

subplot(1,2,2);histogram(DoG.ecc(chan_sel),0:0.1:5,...
    'FaceColor','k','FaceAlpha',0.5);
title('ECC of pRF with good DoG fits & bad LIN fits')
ylabel('nChannels'); xlabel('Ecc');
   
sgtitle('MUA')

if SaveFigs
    saveas(f_neg1,fullfile(figfld, ['EPHYS_NEG-PRF1_MUA.png']));
end
if CloseFigs; close(f_neg1); end
    
%% ---

f_neg3 = figure;
set(f_neg3,'Position',[10 10 1300 1100]);
chan_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh;

subplot(2,2,1); hold on;
plot([0 15],[0 15],'r');
scatter(DoG.ecc(chan_sel),lin.ecc(chan_sel),...
    'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
xlabel('Ecc. DOG')
ylabel('Ecc. LINEAR POS')
set(gca,'xlim',[0 15],'ylim',[0 15]);
yy=get(gca,'ylim');
title('Eccentricity')

subplot(2,2,2); hold on;
histogram(DoG.normamp(chan_sel),-2:.1:2,'FaceColor','k','FaceAlpha',0.5);
xlabel('INH nAMP');ylabel('nChannels');
set(gca,'xlim',[-2 2]);
MM=median(DoG.normamp(chan_sel));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+30]);
title('NORMAMP')

fprintf(['Median nAMP: ' num2str(median(DoG.normamp(chan_sel))) ', IQR: '...
    num2str(iqr(DoG.normamp(chan_sel))) '\n'])

% stats nAmp > 0
% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(DoG.normamp(chan_sel),0,'tail','right');
fprintf(['nAmp > 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

subplot(2,2,3); hold on;
bb = [DoG.ecc(chan_sel) lin.ecc(chan_sel)];
bb2 = [DoG.ang(chan_sel) lin.ang(chan_sel)];
plot([1 2],bb)
plot([1 2],mean(bb),'k','Linewidth',5)
set(gca,'xtick',1:2,'xticklabels',{'DoG','LIN'},...
    'ylim',[0 20],'xlim',[0.8 2.2])
ylabel('Eccentricity');
title('Ecc Diff')

% stats ecc
% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(bb(:,1),bb(:,2));
fprintf(['dEcc ~= 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['\nMean dEcc ' num2str(mean(bb(:,1)-bb(:,2))) ' +/- std ' num2str(std(bb(:,1)-bb(:,2)))])
fprintf(['\nMedian dEcc ' num2str(median(bb(:,1)-bb(:,2))) '\n'])

% stats ang
[p,h,stats] = signrank(bb2(:,1),bb2(:,2));
fprintf(['dAng ~= 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['\nMean dAng ' num2str(mean(bb2(:,1)-bb2(:,2))) ' +/- std ' num2str(std(bb2(:,1)-bb2(:,2)))])
fprintf(['\nMedian dAng ' num2str(median(bb2(:,1)-bb2(:,2))) '\n'])

% mean position diff
dp = sqrt(...
    (DoG.X(chan_sel)-lin.X(chan_sel)).^2 + ...
    (DoG.Y(chan_sel)-lin.Y(chan_sel)).^2);

fprintf(['\nMean dPOS ' num2str(mean(dp)) ' +/- std ' num2str(std(dp)) '\n'])
fprintf(['\nMedian dPOS ' num2str(median(dp)) ' IQR ' num2str(iqr(dp)) '\n'])

subplot(2,2,4); hold on;
histogram(lin.ecc(chan_sel)-DoG.ecc(chan_sel),-20:0.5:20,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Ecc. Diff (POS-DoG)');ylabel('nChannels');
set(gca,'xlim',[-10 10]);
MM=median(bb(:,2)-bb(:,1));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10]);
title('Ecc Diff')

sgtitle(['MUA: pRFs POS LINEAR vs DoG model']);

if SaveFigs
    saveas(f_neg3,fullfile(figfld,['EPHYS_NEG-PRF3_MUA.png']));
end
%if CloseFigs; close(f_neg3); end
    
%% What's specific about the good DoG fits LFP EDITION ====================
% this doesn't happen for MUA
% it is present in LFP (lower freq)

% - find these channels
% - plot their location
% - plot their size
% - plot their prf profile

R2th = 20; % minimum R2
R2enh = 5; % R2 improvement

fb = {'Alpha','Beta'};

for fidx = 1:length(fb)
    
    DoG = tLFP(...
        strcmp(tLFP.Model,'dog_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);
    lin_n = tLFP(...
        strcmp(tLFP.Model,'linear_ephys_cv1_neggain') & strcmp(tLFP.SigType,fb{fidx}),:);
    lin = tLFP(...
        strcmp(tLFP.Model,'linear_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);
    
    %% ---
    f_neg1 = figure;
    set(f_neg1,'Position',[10 10 1200 1000]);
    chan_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh;    
    subplot(2,2,1);scatter(DoG.X(chan_sel),DoG.Y(chan_sel),...
        'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
    set(gca,'xaxislocation','origin','yaxislocation','origin',...
        'xlim',[-8 8],'ylim',[-8 8]);
    title('Locations of pRF with good DoG fits & bad LIN fits')
    xlabel('X deg');ylabel('Y deg');
    
    subplot(2,2,2);histogram(DoG.ecc(chan_sel),0:0.1:5,...
        'FaceColor','k','FaceAlpha',0.5);
    title('ECC of pRF with good DoG fits & bad LIN fits')
    ylabel('nChannels'); xlabel('Ecc');
    
    chan_sel = lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh;
    subplot(2,2,3);scatter(lin_n.X(chan_sel),lin_n.Y(chan_sel),...
        'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
    set(gca,'xaxislocation','origin','yaxislocation','origin',...
        'xlim',[-8 8],'ylim',[-8 8]);
    title('Locations of pRF with good LIN-N fits & bad LIN fits')
    xlabel('X deg');ylabel('Y deg');
    
    subplot(2,2,4);histogram(lin_n.ecc(chan_sel),0:0.1:5,...
        'FaceColor','k','FaceAlpha',0.5);
    title('ECC of pRF with good LIN-N fits & bad LIN fits')
    ylabel('nChannelss'); xlabel('Ecc');
    
    sgtitle(fb{fidx})
    
    if SaveFigs
        saveas(f_neg1,fullfile(figfld, ['EPHYS_NEG-PRF1_' fb{fidx} '.png']));
    end
    if CloseFigs; close(f_neg1); end
    
    %% ---
    f_neg2 = figure;
    set(f_neg2,'Position',[10 10 1300 1600]);
    chan_sel = lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh;
    
    subplot(3,2,1); hold on;
    plot([0 15],[0 15],'r');
    scatter(lin_n.ecc(chan_sel),lin.ecc(chan_sel),...
        'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
    xlabel('Ecc. LINEAR POSNEG')
    ylabel('Ecc. LINEAR POS')
    set(gca,'xlim',[0 15],'ylim',[0 15]);
    yy=get(gca,'ylim');
    title('Eccentricity')
    
    subplot(3,2,2); hold on;
    histogram(lin_n.gain(chan_sel),-1000:25:1000,'FaceColor','k','FaceAlpha',0.5);
    xlabel('gain LIN-POSNEG');ylabel('nChannels');
    set(gca,'xlim',[-1000 1000]);
    MM=median(lin_n.gain(chan_sel));
    yy=get(gca,'ylim');
    plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
    set(gca,'ylim',[0 yy(2)+10]);
    title('Gain')
    
    subplot(3,2,3); hold on;
    bb = [lin_n.ecc(chan_sel) lin.ecc(chan_sel)];
    plot([1 2],bb)
    plot([1 2],mean(bb),'k','Linewidth',5)
    set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
        'ylim',[0 20],'xlim',[0.8 2.2])
    ylabel('Eccentricity');
    title('Ecc Diff')
    
    subplot(3,2,4); hold on;
    histogram(lin.ecc(chan_sel)-lin_n.ecc(chan_sel),-10:0.5:20,...
        'FaceColor','k','FaceAlpha',0.5);
    xlabel('Ecc. Diff (POS-POSNEG)');ylabel('nChannels');
    set(gca,'xlim',[-5 20]);
    MM=median(bb(:,2)-bb(:,1));
    yy=get(gca,'ylim');
    plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
    set(gca,'ylim',[0 yy(2)+10]);
    title('Ecc Diff')
    
    subplot(3,2,5); hold on;
    bb = [lin_n.rfs(chan_sel) lin.rfs(chan_sel)];
    plot([1 2],bb)
    plot([1 2],mean(bb),'k','Linewidth',5)
    set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
        'ylim',[0 6],'xlim',[0.8 2.2])
    ylabel('Size');
    title('Size Diff')
    
    subplot(3,2,6); hold on;
    histogram(lin.rfs(chan_sel)-lin_n.rfs(chan_sel),-10:0.5:10,...
        'FaceColor','k','FaceAlpha',0.5);
    xlabel('Size Diff (POS-POSNEG)');ylabel('nChannels');
    MM=median(bb(:,2)-bb(:,1));
    yy=get(gca,'ylim');
    plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
    set(gca,'ylim',[0 yy(2)+10]);
    set(gca,'xlim',[-5 5]);
    title('Size Diff')
    
    sgtitle([ fb{fidx} ': pRFs POS LINEAR vs POSNEG LINEAR model']);
    
    if SaveFigs 
        saveas(f_neg2,fullfile(figfld,['EPHYS_NEG-PRF2_' fb{fidx} '.png']));
    end
    if CloseFigs; close(f_neg2); end
    
    %% ----
    f_neg3 = figure;
    set(f_neg3,'Position',[10 10 1300 1100]);
    chan_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh & ...
        DoG.normamp~=0 & DoG.ecc<16;
    
    subplot(2,2,1); hold on;
    plot([0 15],[0 15],'r');
    scatter(DoG.ecc(chan_sel),lin.ecc(chan_sel),...
        'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.5);
    xlabel('Ecc. DOG')
    ylabel('Ecc. LINEAR POS')
    set(gca,'xlim',[0 15],'ylim',[0 15]);
    yy=get(gca,'ylim');
    title('Eccentricity')
    
    subplot(2,2,2); hold on;
    histogram(DoG.normamp(chan_sel),-10:5:100,'FaceColor','k','FaceAlpha',0.5);
    xlabel('INH nAMP');ylabel('nChannels');
    set(gca,'xlim',[-10 100]);
    MM=median(DoG.normamp(chan_sel));
    yy=get(gca,'ylim');
    plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
    set(gca,'ylim',[0 yy(2)+30]);
    title('NORMAMP')
    
    subplot(2,2,3); hold on;
    bb = [DoG.ecc(chan_sel) lin.ecc(chan_sel)];
    plot([1 2],bb)
    plot([1 2],mean(bb),'k','Linewidth',5)
    set(gca,'xtick',1:2,'xticklabels',{'DoG','LIN'},...
        'ylim',[0 20],'xlim',[0.8 2.2])
    ylabel('Eccentricity');
    title('Ecc Diff')
    
    subplot(2,2,4); hold on;
    histogram(lin.ecc(chan_sel)-DoG.ecc(chan_sel),-20:0.5:20,...
        'FaceColor','k','FaceAlpha',0.5);
    xlabel('Ecc. Diff (POS-DoG)');ylabel('nChannels');
    set(gca,'xlim',[-5 20]);
    MM=median(bb(:,2)-bb(:,1));
    yy=get(gca,'ylim');
    plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
    set(gca,'ylim',[0 yy(2)+10]);
    title('Ecc Diff')
    
    sgtitle([ fb{fidx} ': pRFs POS LINEAR vs DoG model']);
    
    if SaveFigs
        saveas(f_neg3,fullfile(figfld,['EPHYS_NEG-PRF3_' fb{fidx} '.png']));
    end
    if CloseFigs; close(f_neg3); end
    
end

%% manuscript version ----

R2th = 20; % minimum R2
R2enh = 5; % R2 improvement

fb = {'Alpha','Beta'};

%% ALPHA
fidx = 1;

DoG = tLFP(...
    strcmp(tLFP.Model,'dog_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);
lin_n = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1_neggain') & strcmp(tLFP.SigType,fb{fidx}),:);
lin = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);
css = tLFP(...
    strcmp(tLFP.Model,'css_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);

lin_nGAM1 = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1_neggain') & strcmp(tLFP.SigType,'lGamma'),:);
lin_nGAM2 = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1_neggain') & strcmp(tLFP.SigType,'hGamma'),:);


%% U-LIN ----
f_neg2 = figure; 
roi=1; % only do this for V1 channels
set(f_neg2,'Position',[10 10 900 1200]);
chan_sel = lin_n.Area==roi & lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh;
chan_sel2 = lin_n.Area==roi & lin_n.R2>R2th;

chan_pgain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain>0;
chan_ngain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain<0;

chan_gam1 = lin_n.Area==roi & lin_nGAM1.R2>R2th;
chan_gam2 = lin_n.Area==roi & lin_nGAM2.R2>R2th;

figure(f_neg2);
subplot(3,2,1); hold on;
% gain alpha U-LIN
histogram(lin_n.gain(chan_sel2),-2000:50:2000,'FaceColor','k','FaceAlpha',0.5);
histogram(lin_n.gain(chan_ngain),-2000:50:2000,'FaceColor','r','FaceAlpha',0.5);
histogram(lin_n.gain(chan_pgain),-2000:50:2000,'FaceColor','b','FaceAlpha',0.5);

xlabel('gain U-LIN - ALL ELEC');ylabel('nChannels');
set(gca,'xlim',[-800 1800],'TickDir','out');
MM=median(lin_n.gain(chan_sel2));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10]);
title('Gain')

fprintf(['UNSELECTED - ALPHA MEDIAN GAIN: ' num2str(MM) ', IQR ' num2str(iqr(lin_n.gain(chan_sel))) '\n'])

subplot(3,2,2); hold on;
% gain alpha U-LIN
histogram(lin_n.gain(chan_sel),-1000:50:1000,'FaceColor','k','FaceAlpha',0.5);
xlabel('gain LIN-POSNEG');ylabel('nChannels');
set(gca,'xlim',[-800 1700],'TickDir','out');
MM=median(lin_n.gain(chan_sel));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10]);
title('Gain selected channels (based on R2 ULIN>PLIN')

fprintf(['ALPHA MEDIAN GAIN: ' num2str(MM) ', IQR ' num2str(iqr(lin_n.gain(chan_sel))) '\n'])

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(lin_n.gain(chan_sel),0,'tail','left');
fprintf(['Gain < 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

subplot(3,2,3); hold on;
bb = [lin_n.ecc(chan_sel) lin.ecc(chan_sel)];

mEcc = [mean(lin_n.ecc(chan_ngain)) mean(lin.ecc(chan_ngain)) ...
    mean(lin_n.ecc(chan_pgain))  mean(lin.ecc(chan_pgain)) ];
sdEcc = [std(lin_n.ecc(chan_ngain)) std(lin.ecc(chan_ngain)) ...
    std(lin_n.ecc(chan_pgain))  std(lin.ecc(chan_pgain)) ];

mdEcc = [median(lin_n.ecc(chan_ngain)) median(lin.ecc(chan_ngain)) ...
    median(lin_n.ecc(chan_pgain))  median(lin.ecc(chan_pgain)) ];
iqrEcc = [iqr(lin_n.ecc(chan_ngain)) iqr(lin.ecc(chan_ngain)) ...
    iqr(lin_n.ecc(chan_pgain))  iqr(lin.ecc(chan_pgain)) ]./2;

% plot([1 2],bb)
% plot([1 2],mean(bb),'k','Linewidth',5)
% errorbar([1 2],mean(bb),std(bb),...
%     'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
% set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
%     'ylim',[0 30],'xlim',[0.8 2.2],'TickDir','out')

% errorbar(1:4,mEcc,sdEcc,...
%     'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
errorbar(1:4,mdEcc,iqrEcc,...
    'ko','MarkerSize',4,'MarkerFaceColor','k','Linewidth',2)
set(gca,'xtick',1:4,'xticklabels',{'ULIN-N','PLIN-N','ULIN-P','PLIN-P'},...
    'ylim',[0 25],'xlim',[0.8 4.2],'TickDir','out')

ylabel('Eccentricity');
title('Ecc')

% Ecc difference
fprintf('-- STATS Ecc --\n');
fprintf(['mEcc nGain U-LIN: ' num2str(mEcc(1)) ', STD ' num2str(sdEcc(1)) '\n'])
fprintf(['mEcc pGain U-LIN: ' num2str(mEcc(3)) ', STD ' num2str(sdEcc(3)) '\n'])
fprintf(['mdEcc nGain U-LIN: ' num2str(mdEcc(1)) ', IQR ' num2str(iqrEcc(1)) '\n'])
fprintf(['mdEcc pGain U-LIN: ' num2str(mdEcc(3)) ', IQR ' num2str(iqrEcc(3)) '\n'])
[p,h,stats] = ranksum(lin_n.ecc(chan_ngain),lin_n.ecc(chan_pgain));
fprintf(['Ecc difference pos gain U-LIN vs neg gain U-LIN: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

fprintf(['mEcc nGain U-LIN: ' num2str(mEcc(1)) ', STD ' num2str(sdEcc(1)) '\n'])
fprintf(['mEcc nGain P-LIN: ' num2str(mEcc(2)) ', STD ' num2str(sdEcc(2)) '\n'])
fprintf(['mdEcc nGain U-LIN: ' num2str(mdEcc(1)) ', STD ' num2str(iqrEcc(1)) '\n'])
fprintf(['mdEcc nGain P-LIN: ' num2str(mdEcc(2)) ', STD ' num2str(iqrEcc(2)) '\n'])
[p,h,stats] = signrank(lin_n.ecc(chan_ngain),lin.ecc(chan_ngain));
fprintf(['Ecc difference neg gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

fprintf(['mEcc pGain U-LIN: ' num2str(mEcc(3)) ', STD ' num2str(sdEcc(3)) '\n'])
fprintf(['mEcc pGain P-LIN: ' num2str(mEcc(4)) ', STD ' num2str(sdEcc(4)) '\n'])
fprintf(['mdEcc pGain U-LIN: ' num2str(mdEcc(3)) ', STD ' num2str(iqrEcc(3)) '\n'])
fprintf(['mdEcc pGain P-LIN: ' num2str(mdEcc(4)) ', STD ' num2str(iqrEcc(4)) '\n'])
[p,h,stats] = ranksum(lin_n.ecc(chan_pgain),lin.ecc(chan_pgain));
fprintf(['Ecc difference pos gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

subplot(3,2,4); hold on;
% histogram(lin.ecc(chan_sel)-lin_n.ecc(chan_sel),-50:1:50,...
%     'FaceColor','k','FaceAlpha',0.5);
histogram(lin.ecc(chan_ngain)-lin_n.ecc(chan_ngain),-50:1:50,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Ecc. Diff (PLIN-ULIN_N)','interpreter','none');ylabel('nChannels');
set(gca,'xlim',[-10 35],'TickDir','out');

MM=median(lin.ecc(chan_ngain)-lin_n.ecc(chan_ngain));

yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+5]);
title('Ecc Diff')

% Wilcoxon 
[p,h,stats] = signrank(bb(:,1),bb(:,2));
fprintf(['ECC diff Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

fprintf(['ALPHA MEDIAN ECC DIFF: ' num2str(median(lin.ecc(chan_sel)-lin_n.ecc(chan_sel))) ...
    ', IQR ' num2str(iqr(lin.ecc(chan_sel)-lin_n.ecc(chan_sel))) '\n'])

subplot(3,2,5); hold on;
bb = [lin_n.rfs(chan_sel) lin.rfs(chan_sel)];

mSz = [mean(lin_n.rfs(chan_ngain)) mean(lin.rfs(chan_ngain)) ...
    mean(lin_n.rfs(chan_pgain))  mean(lin.rfs(chan_pgain)) ];
sdSz = [std(lin_n.rfs(chan_ngain)) std(lin.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_pgain))  std(lin.rfs(chan_pgain)) ];

mdSz = [median(lin_n.rfs(chan_ngain)) median(lin.rfs(chan_ngain)) ...
    median(lin_n.rfs(chan_pgain))  median(lin.rfs(chan_pgain)) ];
iqrSz = [iqr(lin_n.rfs(chan_ngain)) iqr(lin.rfs(chan_ngain)) ...
    iqr(lin_n.rfs(chan_pgain))  iqr(lin.rfs(chan_pgain)) ]./2;

%plot([1 2],bb)
%plot([1 2],mean(bb),'k','Linewidth',5)
% errorbar([1 2],mean(bb),std(bb),'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
% set(gca,'xtick',1:2,'xticklabels',{'U-LIN','P-LIN'},...
%     'ylim',[0 6],'xlim',[0.8 2.2],'TickDir','out')
%errorbar(1:4,mSz,sdSz,'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
errorbar(1:4,mdSz,iqrSz,'ko','MarkerSize',4,'MarkerFaceColor','k','Linewidth',2)
set(gca,'xtick',1:4,'xticklabels',{'ULIN-N','PLIN-N','ULIN-P','PLIN-P'},...
    'ylim',[0 6],'xlim',[0.8 4.2],'TickDir','out')
ylabel('Size');
title('Size')


% Sz difference
fprintf('-- STATS Sz --\n');
fprintf(['mSz nGain U-LIN: ' num2str(mSz(1)) ', STD ' num2str(sdSz(1)) '\n'])
fprintf(['mSz pGain U-LIN: ' num2str(mSz(3)) ', STD ' num2str(sdSz(3)) '\n'])
fprintf(['mdSz nGain U-LIN: ' num2str(mdSz(1)) ', IQR ' num2str(iqrSz(1)) '\n'])
fprintf(['mdSz pGain U-LIN: ' num2str(mdSz(3)) ', IQR ' num2str(iqrSz(3)) '\n'])
[p,h,stats] = ranksum(lin_n.rfs(chan_ngain),lin_n.rfs(chan_pgain));
fprintf(['Sz difference pos gain U-LIN vs neg gain U-LIN: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['mSz nGain U-LIN: ' num2str(mSz(1)) ', STD ' num2str(sdSz(1)) '\n'])
fprintf(['mSz nGain P-LIN: ' num2str(mSz(2)) ', STD ' num2str(sdSz(2)) '\n'])
fprintf(['mdSz nGain U-LIN: ' num2str(mdSz(1)) ', IQR ' num2str(iqrSz(1)) '\n'])
fprintf(['mdSz nGain P-LIN: ' num2str(mdSz(2)) ', IQR ' num2str(iqrSz(2)) '\n'])
[p,h,stats] = ranksum(lin_n.rfs(chan_ngain),lin.rfs(chan_ngain));
fprintf(['Sz difference neg gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['mSz pGain U-LIN: ' num2str(mSz(3)) ', STD ' num2str(sdSz(3)) '\n'])
fprintf(['mSz pGain P-LIN: ' num2str(mSz(4)) ', STD ' num2str(sdSz(4)) '\n'])
fprintf(['mdSz pGain U-LIN: ' num2str(mdSz(3)) ', IQR ' num2str(iqrSz(3)) '\n'])
fprintf(['mdSz pGain P-LIN: ' num2str(mdSz(4)) ', IQR ' num2str(iqrSz(4)) '\n'])
[p,h,stats] = ranksum(lin_n.rfs(chan_pgain),lin.rfs(chan_pgain));
fprintf(['Ecc difference pos gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);



subplot(3,2,6); hold on;
% histogram(lin.rfs(chan_sel)-lin_n.rfs(chan_sel),-10:0.5:10,...
%     'FaceColor','k','FaceAlpha',0.5);
histogram(lin.rfs(chan_ngain)-lin_n.rfs(chan_ngain),-10:0.5:10,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Size Diff (POS-POSNEG)');ylabel('nChannels');
% MM=median(bb(:,2)-bb(:,1));
MM=median(lin.rfs(chan_ngain)-lin_n.rfs(chan_ngain));

yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10]);
set(gca,'xlim',[-5 10],'TickDir','out');
title('Size Diff')

% Wilcoxon 
[p,h,stats] = signrank(bb(:,1),bb(:,2));
fprintf(['SZ diff Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

sgtitle([ fb{fidx} ': pRFs P-LIN vs U-LIN']);


%%
fextra=figure;
subplot(2,2,1);hold on;
plot([0 5],[0 5]);
scatter(lin_n.ecc(chan_ngain),lin_nGAM1.ecc(chan_ngain));
%set(gca,'xlim',[0 5],'ylim',[0 5]);
subplot(2,2,2);hold on;
plot([0 5],[0 5]);
scatter(lin_n.ecc(chan_pgain),lin_nGAM1.ecc(chan_pgain));
%set(gca,'xlim',[0 5],'ylim',[0 5]);
subplot(2,2,3);hold on;
plot([0 5],[0 5]);
scatter(lin_n.ecc(chan_ngain),lin_nGAM2.ecc(chan_ngain));
%set(gca,'xlim',[0 5],'ylim',[0 5]);
subplot(2,2,4);hold on;
plot([0 5],[0 5]);
scatter(lin_n.ecc(chan_pgain),lin_nGAM2.ecc(chan_pgain));
%set(gca,'xlim',[0 5],'ylim',[0 5]);

%%
% distance
lin_n = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1_neggain') & strcmp(tLFP.SigType,'Alpha'),:);

roi=1; % only do this for V1 channels
chan_sel = lin_n.Area==roi & lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh;
chan_sel2 = lin_n.Area==roi & lin_n.R2>R2th;

chan_pgain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain>0;
chan_ngain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain<0;

chan_gam1 = lin_n.Area==roi & lin_nGAM1.R2>R2th;
chan_gam2 = lin_n.Area==roi & lin_nGAM2.R2>R2th;
[mean(sqrt((lin_n.X(chan_ngain)-lin_nGAM1.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2)) ...
    std(sqrt((lin_n.X(chan_ngain)-lin_nGAM1.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2))./sqrt(sum(chan_ngain))]
[mean(sqrt((lin_n.X(chan_pgain)-lin_nGAM1.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM1.Y(chan_pgain)).^2)) ...
    std(sqrt((lin_n.X(chan_pgain)-lin_nGAM1.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM1.Y(chan_pgain)).^2))./sqrt(sum(chan_pgain))]

[mean(sqrt((lin_n.X(chan_ngain)-lin_nGAM2.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2)) ...
    std(sqrt((lin_n.X(chan_ngain)-lin_nGAM2.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM2.Y(chan_ngain)).^2))./sqrt(sum(chan_ngain))]
[mean(sqrt((lin_n.X(chan_pgain)-lin_nGAM2.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM2.Y(chan_pgain)).^2)) ...
    std(sqrt((lin_n.X(chan_pgain)-lin_nGAM2.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM2.Y(chan_pgain)).^2))./sqrt(sum(chan_pgain))]


% size
[mean(lin_n.rfs(chan_ngain)-lin_nGAM1.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_ngain)-lin_nGAM1.rfs(chan_ngain))./sqrt(sum(chan_ngain))]
[mean(lin_n.rfs(chan_pgain)-lin_nGAM1.rfs(chan_pgain)) ...
    std(lin_n.rfs(chan_pgain)-lin_nGAM1.rfs(chan_pgain))./sqrt(sum(chan_pgain))]

% normalized size
[mean(lin_n.rfs(chan_ngain)./lin_nGAM1.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_ngain)./lin_nGAM1.rfs(chan_ngain))./sqrt(sum(chan_ngain))]
[mean(lin_n.rfs(chan_pgain)./lin_nGAM1.rfs(chan_pgain)) ...
    std(lin_n.rfs(chan_pgain)./lin_nGAM1.rfs(chan_pgain))./sqrt(sum(chan_pgain))]

[mean(lin_n.rfs(chan_ngain)./lin_nGAM2.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_ngain)./lin_nGAM2.rfs(chan_ngain))./sqrt(sum(chan_ngain))]
[mean(lin_n.rfs(chan_pgain)./lin_nGAM2.rfs(chan_pgain)) ...
    std(lin_n.rfs(chan_pgain)./lin_nGAM2.rfs(chan_pgain))./sqrt(sum(chan_pgain))]

% distance in relation to size
distn = sqrt((lin_n.X(chan_ngain)-lin_nGAM1.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2);
distp = sqrt((lin_n.X(chan_pgain)-lin_nGAM1.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM1.Y(chan_pgain)).^2);
sz2n=lin_n.rfs(chan_ngain)+lin_nGAM1.rfs(chan_ngain);
sz2p=lin_n.rfs(chan_pgain)+lin_nGAM1.rfs(chan_pgain);
figure;hold on;
scatter(distn,sz2n);
scatter(distp,sz2p);
plot([0 20],[0 20]);
set(gca,'xlim',[0 20],'ylim',[0 20]);
legend({'negative', 'positive'});

[mean(distn./sz2n) std(distn./sz2n)./sqrt(sum(chan_ngain))]
[mean(distp./sz2p) std(distp./sz2p)./sqrt(sum(chan_pgain))]

distn = sqrt((lin_n.X(chan_ngain)-lin_nGAM2.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM2.Y(chan_ngain)).^2);
distp = sqrt((lin_n.X(chan_pgain)-lin_nGAM2.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM2.Y(chan_pgain)).^2);
sz2n=lin_n.rfs(chan_ngain)+lin_nGAM2.rfs(chan_ngain);
sz2p=lin_n.rfs(chan_pgain)+lin_nGAM2.rfs(chan_pgain);
figure;hold on;
scatter(distn,sz2n);
scatter(distp,sz2p);
plot([0 20],[0 20]);
set(gca,'xlim',[0 20],'ylim',[0 20]);
legend({'negative', 'positive'});
[mean(distn./sz2n) std(distn./sz2n)./sqrt(sum(chan_ngain))]
[mean(distp./sz2p) std(distp./sz2p)./sqrt(sum(chan_pgain))]



%% Beta
% distance
lin_n = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1_neggain') & strcmp(tLFP.SigType,'Beta'),:);

roi=1; % only do this for V1 channels
chan_sel = lin_n.Area==roi & lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh;
chan_sel2 = lin_n.Area==roi & lin_n.R2>R2th;

chan_pgain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain>0;
chan_ngain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain<0;

chan_gam1 = lin_n.Area==roi & lin_nGAM1.R2>R2th;
chan_gam2 = lin_n.Area==roi & lin_nGAM2.R2>R2th;

[mean(sqrt((lin_n.X(chan_ngain)-lin_nGAM1.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2)) ...
    std(sqrt((lin_n.X(chan_ngain)-lin_nGAM1.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2))./sqrt(sum(chan_ngain))]
[mean(sqrt((lin_n.X(chan_pgain)-lin_nGAM1.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM1.Y(chan_pgain)).^2)) ...
    std(sqrt((lin_n.X(chan_pgain)-lin_nGAM1.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM1.Y(chan_pgain)).^2))./sqrt(sum(chan_pgain))]

[mean(sqrt((lin_n.X(chan_ngain)-lin_nGAM2.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2)) ...
    std(sqrt((lin_n.X(chan_ngain)-lin_nGAM2.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM2.Y(chan_ngain)).^2))./sqrt(sum(chan_ngain))]
[mean(sqrt((lin_n.X(chan_pgain)-lin_nGAM2.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM2.Y(chan_pgain)).^2)) ...
    std(sqrt((lin_n.X(chan_pgain)-lin_nGAM2.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM2.Y(chan_pgain)).^2))./sqrt(sum(chan_pgain))]

% size
[mean(lin_n.rfs(chan_ngain)-lin_nGAM1.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_ngain)-lin_nGAM1.rfs(chan_ngain))./sqrt(sum(chan_ngain))]
[mean(lin_n.rfs(chan_pgain)-lin_nGAM1.rfs(chan_pgain)) ...
    std(lin_n.rfs(chan_pgain)-lin_nGAM1.rfs(chan_pgain))./sqrt(sum(chan_pgain))]

% normalized size
[mean(lin_n.rfs(chan_ngain)./lin_nGAM1.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_ngain)./lin_nGAM1.rfs(chan_ngain))./sqrt(sum(chan_ngain))]
[mean(lin_n.rfs(chan_pgain)./lin_nGAM1.rfs(chan_pgain)) ...
    std(lin_n.rfs(chan_pgain)./lin_nGAM1.rfs(chan_pgain))./sqrt(sum(chan_pgain))]

[mean(lin_n.rfs(chan_ngain)./lin_nGAM2.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_ngain)./lin_nGAM2.rfs(chan_ngain))./sqrt(sum(chan_ngain))]
[mean(lin_n.rfs(chan_pgain)./lin_nGAM2.rfs(chan_pgain)) ...
    std(lin_n.rfs(chan_pgain)./lin_nGAM2.rfs(chan_pgain))./sqrt(sum(chan_pgain))]

% distance in relation to size
distn = sqrt((lin_n.X(chan_ngain)-lin_nGAM1.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM1.Y(chan_ngain)).^2);
distp = sqrt((lin_n.X(chan_pgain)-lin_nGAM1.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM1.Y(chan_pgain)).^2);
sz2n=lin_n.rfs(chan_ngain)+lin_nGAM1.rfs(chan_ngain);
sz2p=lin_n.rfs(chan_pgain)+lin_nGAM1.rfs(chan_pgain);
figure;hold on;
scatter(distn,sz2n);
scatter(distp,sz2p);
plot([0 20],[0 20]);
set(gca,'xlim',[0 20],'ylim',[0 20]);
legend({'negative', 'positive'});

figure;hold on;
[mean(distn./sz2n) std(distn./sz2n)./sqrt(sum(chan_ngain))...
    median(distn./sz2n) iqr(distn./sz2n)./2]
histogram(distn./sz2n,100,'Normalization','probability')
[mean(distp./sz2p) std(distp./sz2p)./sqrt(sum(chan_pgain))...
    median(distp./sz2p) iqr(distp./sz2p)./2]
histogram(distp./sz2p,100,'Normalization','probability')
title('Low Gamma'); 


distn = sqrt((lin_n.X(chan_ngain)-lin_nGAM2.X(chan_ngain)).^2 + ...
    (lin_n.Y(chan_ngain)-lin_nGAM2.Y(chan_ngain)).^2);
distp = sqrt((lin_n.X(chan_pgain)-lin_nGAM2.X(chan_pgain)).^2 + ...
    (lin_n.Y(chan_pgain)-lin_nGAM2.Y(chan_pgain)).^2);
sz2n=lin_n.rfs(chan_ngain)+lin_nGAM2.rfs(chan_ngain);
sz2p=lin_n.rfs(chan_pgain)+lin_nGAM2.rfs(chan_pgain);
figure;hold on;
scatter(distn,sz2n);
scatter(distp,sz2p);
plot([0 20],[0 20]);
set(gca,'xlim',[0 20],'ylim',[0 20]);
legend({'negative', 'positive'});

figure;hold on;
[mean(distn./sz2n) std(distn./sz2n)./sqrt(sum(chan_ngain))...
    median(distn./sz2n) iqr(distn./sz2n)./2]
histogram(distn./sz2n,100,'Normalization','probability')
[mean(distp./sz2p) std(distp./sz2p)./sqrt(sum(chan_pgain))...
    median(distp./sz2p) iqr(distp./sz2p)./2]
histogram(distp./sz2p,100,'Normalization','probability')
title('High Gamma'); 










%% ----

f_neg3 = figure;
set(f_neg3,'Position',[10 10 900 800]);
chan_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh & ...
    DoG.normamp~=0 & DoG.ecc<16;

subplot(2,2,2); hold on;
histogram(DoG.normamp(chan_sel),-10:5:200,'FaceColor','k','FaceAlpha',0.5);
xlabel('INH nAMP');ylabel('nChannels');
set(gca,'xlim',[-10 120]);
MM=median(DoG.normamp(chan_sel));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+30],'TickDir','out');
title('NORMAMP')

fprintf(['ALPHA MEDIAN NAMP: ' num2str(MM) ', IQR ' num2str(iqr(DoG.normamp(chan_sel))) '\n'])

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(lin_n.gain(chan_sel),0,'tail','left');
fprintf(['Gain < 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);


subplot(2,2,3); hold on;
bb = [DoG.ecc(chan_sel) lin.ecc(chan_sel)];
%plot([1 2],bb)
%plot([1 2],mean(bb),'k','Linewidth',5)
errorbar([1 2],mean(bb),std(bb),'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
set(gca,'xtick',1:2,'xticklabels',{'DoG','LIN'},...
    'ylim',[0 25],'xlim',[0.8 2.2],'TickDir','out')
ylabel('Eccentricity');
title('Ecc Diff')

% Wilcoxon 
[p,h,stats] = signrank(bb(:,1),bb(:,2));
fprintf(['ECC diff Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

subplot(2,2,4); hold on;
histogram(lin.ecc(chan_sel)-DoG.ecc(chan_sel),-20:1:30,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Ecc. Diff (POS-DoG)');ylabel('nChannels');
set(gca,'xlim',[-5 30]);
MM=median(bb(:,2)-bb(:,1));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10],'TickDir','out');
title('Ecc Diff')

sgtitle([ fb{fidx} ': pRFs POS LINEAR vs DoG model']);

fprintf(['ALPHA MEDIAN ECC DIFF: ' num2str(median(lin.ecc(chan_sel)-DoG.ecc(chan_sel))) ...
    ', IQR ' num2str(iqr(lin.ecc(chan_sel)-DoG.ecc(chan_sel))) '\n'])


if SaveFigs
    saveas(f_neg3,fullfile(figfld,['EPHYS_NEG-PRF3_' fb{fidx} '.png']));
end
if CloseFigs; close(f_neg3); end
    
%% BETA
fidx = 2;

DoG = tLFP(...
    strcmp(tLFP.Model,'dog_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);
lin_n = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1_neggain') & strcmp(tLFP.SigType,fb{fidx}),:);
lin = tLFP(...
    strcmp(tLFP.Model,'linear_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);
css = tLFP(...
    strcmp(tLFP.Model,'css_ephys_cv1') & strcmp(tLFP.SigType,fb{fidx}),:);

%% ----
f_neg2 = figure;
roi=1; % only do this for V1 channels
set(f_neg2,'Position',[10 10 900 1200]);
chan_sel = lin_n.Area==roi & lin_n.R2>R2th & lin_n.R2>lin.R2+R2enh;
chan_sel2 = lin_n.Area==roi & lin_n.R2>R2th;

chan_pgain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain>0;
chan_ngain = lin_n.Area==roi & lin_n.R2>R2th & lin_n.gain<0;

subplot(3,2,1); hold on;
% gain alpha U-LIN
histogram(lin_n.gain(chan_sel2),-2000:50:2000,'FaceColor','k','FaceAlpha',0.5);
histogram(lin_n.gain(chan_ngain),-2000:50:2000,'FaceColor','r','FaceAlpha',0.5);
histogram(lin_n.gain(chan_pgain),-2000:50:2000,'FaceColor','b','FaceAlpha',0.5);

xlabel('gain U-LIN - ALL ELEC');ylabel('nChannels');
set(gca,'xlim',[-800 1800],'TickDir','out');
MM=median(lin_n.gain(chan_sel2));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10]);
title('Gain')

fprintf(['UNSELECTED - BETA MEDIAN GAIN: ' num2str(MM) ', IQR ' num2str(iqr(lin_n.gain(chan_sel))) '\n'])

subplot(3,2,2); hold on;
% gain alpha U-LIN
histogram(lin_n.gain(chan_sel),-1000:50:1000,'FaceColor','k','FaceAlpha',0.5);
xlabel('gain LIN-POSNEG');ylabel('nChannels');
set(gca,'xlim',[-800 1700],'TickDir','out');
MM=median(lin_n.gain(chan_sel));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10]);
title('Gain selected channels (based on R2 ULIN>PLIN')

fprintf(['BETA MEDIAN GAIN: ' num2str(MM) ', IQR ' num2str(iqr(lin_n.gain(chan_sel))) '\n'])

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(lin_n.gain(chan_sel),0,'tail','left');
fprintf(['Gain < 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

subplot(3,2,3); hold on;
bb = [lin_n.ecc(chan_sel) lin.ecc(chan_sel)];

mEcc = [mean(lin_n.ecc(chan_ngain)) mean(lin.ecc(chan_ngain)) ...
    mean(lin_n.ecc(chan_pgain))  mean(lin.ecc(chan_pgain)) ];
sdEcc = [std(lin_n.ecc(chan_ngain)) std(lin.ecc(chan_ngain)) ...
    std(lin_n.ecc(chan_pgain))  std(lin.ecc(chan_pgain)) ];

mdEcc = [median(lin_n.ecc(chan_ngain)) median(lin.ecc(chan_ngain)) ...
    median(lin_n.ecc(chan_pgain))  median(lin.ecc(chan_pgain)) ];
iqrEcc = [iqr(lin_n.ecc(chan_ngain)) iqr(lin.ecc(chan_ngain)) ...
    iqr(lin_n.ecc(chan_pgain))  iqr(lin.ecc(chan_pgain)) ]./2;

% plot([1 2],bb)
% plot([1 2],mean(bb),'k','Linewidth',5)
% errorbar([1 2],mean(bb),std(bb),...
%     'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
% set(gca,'xtick',1:2,'xticklabels',{'LIN-N','LIN'},...
%     'ylim',[0 30],'xlim',[0.8 2.2],'TickDir','out')

% errorbar(1:4,mEcc,sdEcc,...
%     'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
errorbar(1:4,mdEcc,iqrEcc,...
    'ko','MarkerSize',4,'MarkerFaceColor','k','Linewidth',2)
set(gca,'xtick',1:4,'xticklabels',{'ULIN-N','PLIN-N','ULIN-P','PLIN-P'},...
    'ylim',[0 35],'xlim',[0.8 4.2],'TickDir','out')

ylabel('Eccentricity');
title('Ecc')

% Ecc difference
fprintf('-- STATS Ecc --\n');
fprintf(['mEcc nGain U-LIN: ' num2str(mEcc(1)) ', STD ' num2str(sdEcc(1)) '\n'])
fprintf(['mEcc pGain U-LIN: ' num2str(mEcc(3)) ', STD ' num2str(sdEcc(3)) '\n'])
fprintf(['mdEcc nGain U-LIN: ' num2str(mdEcc(1)) ', IQR ' num2str(iqrEcc(1)) '\n'])
fprintf(['mdEcc pGain U-LIN: ' num2str(mdEcc(3)) ', IQR ' num2str(iqrEcc(3)) '\n'])
[p,h,stats] = ranksum(lin_n.ecc(chan_ngain),lin_n.ecc(chan_pgain));
fprintf(['Ecc difference pos gain U-LIN vs neg gain U-LIN: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['mEcc nGain U-LIN: ' num2str(mEcc(1)) ', STD ' num2str(sdEcc(1)) '\n'])
fprintf(['mEcc nGain P-LIN: ' num2str(mEcc(2)) ', STD ' num2str(sdEcc(2)) '\n'])
fprintf(['mdEcc nGain U-LIN: ' num2str(mdEcc(1)) ', STD ' num2str(iqrEcc(1)) '\n'])
fprintf(['mdEcc nGain P-LIN: ' num2str(mdEcc(2)) ', STD ' num2str(iqrEcc(2)) '\n'])
[p,h,stats] = ranksum(lin_n.ecc(chan_ngain),lin.ecc(chan_ngain));
fprintf(['Ecc difference neg gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['mEcc pGain U-LIN: ' num2str(mEcc(3)) ', STD ' num2str(sdEcc(3)) '\n'])
fprintf(['mEcc pGain P-LIN: ' num2str(mEcc(4)) ', STD ' num2str(sdEcc(4)) '\n'])
fprintf(['mdEcc pGain U-LIN: ' num2str(mdEcc(3)) ', STD ' num2str(iqrEcc(3)) '\n'])
fprintf(['mdEcc pGain P-LIN: ' num2str(mdEcc(4)) ', STD ' num2str(iqrEcc(4)) '\n'])
[p,h,stats] = ranksum(lin_n.ecc(chan_pgain),lin.ecc(chan_pgain));
fprintf(['Ecc difference pos gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

subplot(3,2,4); hold on;
% histogram(lin.ecc(chan_sel)-lin_n.ecc(chan_sel),-50:1:50,...
%     'FaceColor','k','FaceAlpha',0.5);
histogram(lin.ecc(chan_ngain)-lin_n.ecc(chan_ngain),-50:1:50,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Ecc. Diff (PLIN-ULIN_N)','interpreter','none');ylabel('nChannels');
set(gca,'xlim',[-10 35],'TickDir','out');

MM=median(lin.ecc(chan_ngain)-lin_n.ecc(chan_ngain));

yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+5]);
title('Ecc Diff')

% Wilcoxon 
[p,h,stats] = signrank(bb(:,1),bb(:,2));
fprintf(['ECC diff Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

fprintf(['BETA MEDIAN ECC DIFF: ' num2str(median(lin.ecc(chan_sel)-lin_n.ecc(chan_sel))) ...
    ', IQR ' num2str(iqr(lin.ecc(chan_sel)-lin_n.ecc(chan_sel))) '\n'])

subplot(3,2,5); hold on;
bb = [lin_n.rfs(chan_sel) lin.rfs(chan_sel)];

mSz = [mean(lin_n.rfs(chan_ngain)) mean(lin.rfs(chan_ngain)) ...
    mean(lin_n.rfs(chan_pgain))  mean(lin.rfs(chan_pgain)) ];
sdSz = [std(lin_n.rfs(chan_ngain)) std(lin.rfs(chan_ngain)) ...
    std(lin_n.rfs(chan_pgain))  std(lin.rfs(chan_pgain)) ];

mdSz = [median(lin_n.rfs(chan_ngain)) median(lin.rfs(chan_ngain)) ...
    median(lin_n.rfs(chan_pgain))  median(lin.rfs(chan_pgain)) ];
iqrSz = [iqr(lin_n.rfs(chan_ngain)) iqr(lin.rfs(chan_ngain)) ...
    iqr(lin_n.rfs(chan_pgain))  iqr(lin.rfs(chan_pgain)) ]./2;

%plot([1 2],bb)
%plot([1 2],mean(bb),'k','Linewidth',5)
% errorbar([1 2],mean(bb),std(bb),'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
% set(gca,'xtick',1:2,'xticklabels',{'U-LIN','P-LIN'},...
%     'ylim',[0 6],'xlim',[0.8 2.2],'TickDir','out')
%errorbar(1:4,mSz,sdSz,'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
errorbar(1:4,mdSz,iqrSz,'ko','MarkerSize',4,'MarkerFaceColor','k','Linewidth',2)
set(gca,'xtick',1:4,'xticklabels',{'ULIN-N','PLIN-N','ULIN-P','PLIN-P'},...
    'ylim',[0 6],'xlim',[0.8 4.2],'TickDir','out')
ylabel('Size');
title('Size')


% Sz difference
fprintf('-- STATS Sz --\n');
fprintf(['mSz nGain U-LIN: ' num2str(mSz(1)) ', STD ' num2str(sdSz(1)) '\n'])
fprintf(['mSz pGain U-LIN: ' num2str(mSz(3)) ', STD ' num2str(sdSz(3)) '\n'])
fprintf(['mdSz nGain U-LIN: ' num2str(mdSz(1)) ', IQR ' num2str(iqrSz(1)) '\n'])
fprintf(['mdSz pGain U-LIN: ' num2str(mdSz(3)) ', IQR ' num2str(iqrSz(3)) '\n'])
[p,h,stats] = ranksum(lin_n.rfs(chan_ngain),lin_n.rfs(chan_pgain));
fprintf(['Sz difference pos gain U-LIN vs neg gain U-LIN: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['mSz nGain U-LIN: ' num2str(mSz(1)) ', STD ' num2str(sdSz(1)) '\n'])
fprintf(['mSz nGain P-LIN: ' num2str(mSz(2)) ', STD ' num2str(sdSz(2)) '\n'])
fprintf(['mdSz nGain U-LIN: ' num2str(mdSz(1)) ', IQR ' num2str(iqrSz(1)) '\n'])
fprintf(['mdSz nGain P-LIN: ' num2str(mdSz(2)) ', IQR ' num2str(iqrSz(2)) '\n'])
[p,h,stats] = ranksum(lin_n.rfs(chan_ngain),lin.rfs(chan_ngain));
fprintf(['Sz difference neg gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);
fprintf(['mSz pGain U-LIN: ' num2str(mSz(3)) ', STD ' num2str(sdSz(3)) '\n'])
fprintf(['mSz pGain P-LIN: ' num2str(mSz(4)) ', STD ' num2str(sdSz(4)) '\n'])
fprintf(['mdSz pGain U-LIN: ' num2str(mdSz(3)) ', IQR ' num2str(iqrSz(3)) '\n'])
fprintf(['mdSz pGain P-LIN: ' num2str(mdSz(4)) ', IQR ' num2str(iqrSz(4)) '\n'])
[p,h,stats] = ranksum(lin_n.rfs(chan_pgain),lin.rfs(chan_pgain));
fprintf(['Sz difference pos gains U-LIN vs P-LIN : Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);



subplot(3,2,6); hold on;
% histogram(lin.rfs(chan_sel)-lin_n.rfs(chan_sel),-10:0.5:10,...
%     'FaceColor','k','FaceAlpha',0.5);
histogram(lin.rfs(chan_ngain)-lin_n.rfs(chan_ngain),-10:0.5:10,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Size Diff (POS-POSNEG)');ylabel('nChannels');
% MM=median(bb(:,2)-bb(:,1));
MM=median(lin.rfs(chan_ngain)-lin_n.rfs(chan_ngain));

yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10]);
set(gca,'xlim',[-5 10],'TickDir','out');
title('Size Diff')

% Wilcoxon 
[p,h,stats] = signrank(bb(:,1),bb(:,2));
fprintf(['SZ diff Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

sgtitle([ fb{fidx} ': pRFs P-LIN vs U-LIN']);

%% ----

f_neg3 = figure;
set(f_neg3,'Position',[10 10 900 800],'Renderer','painters');
chan_sel = DoG.R2>R2th & DoG.R2>lin.R2+R2enh & ...
    DoG.normamp~=0 & DoG.ecc<16;

subplot(2,2,2); hold on;
histogram(DoG.normamp(chan_sel),-10:2:50,'FaceColor','k','FaceAlpha',0.5);
xlabel('INH nAMP');ylabel('nChannels');
set(gca,'xlim',[-5 35]);
MM=median(DoG.normamp(chan_sel));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+30],'TickDir','out');
title('NORMAMP')

fprintf(['BETA MEDIAN NAMP: ' num2str(MM) ', IQR ' num2str(iqr(DoG.normamp(chan_sel))) '\n'])

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(lin_n.gain(chan_sel),0,'tail','left');
fprintf(['Gain < 0: Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);


subplot(2,2,3); hold on;
bb = [DoG.ecc(chan_sel) lin.ecc(chan_sel)];
%plot([1 2],bb)
%plot([1 2],mean(bb),'k','Linewidth',5)
errorbar([1 2],mean(bb),std(bb),'ko','MarkerSize',10,'MarkerFaceColor','k','Linewidth',2)
set(gca,'xtick',1:2,'xticklabels',{'DoG','LIN'},...
    'ylim',[0 30],'xlim',[0.8 2.2],'TickDir','out')
ylabel('Eccentricity');
title('Ecc Diff')

% Wilcoxon 
[p,h,stats] = signrank(bb(:,1),bb(:,2));
fprintf(['ECC diff Wilcoxon z = ' ...
    num2str(stats.zval) ', p = ' num2str(p) '\n']);

subplot(2,2,4); hold on;
histogram(lin.ecc(chan_sel)-DoG.ecc(chan_sel),-20:1:30,...
    'FaceColor','k','FaceAlpha',0.5);
xlabel('Ecc. Diff (POS-DoG)');ylabel('nChannels');
set(gca,'xlim',[-18 30]);
MM=median(bb(:,2)-bb(:,1));
yy=get(gca,'ylim');
plot([MM MM], [0 yy(2)+40],'k','Linewidth',5)
set(gca,'ylim',[0 yy(2)+10],'TickDir','out');
title('Ecc Diff')

sgtitle([ fb{fidx} ': pRFs POS LINEAR vs DoG model']);

fprintf(['BETA MEDIAN ECC DIFF: ' num2str(median(lin.ecc(chan_sel)-DoG.ecc(chan_sel))) ...
    ', IQR ' num2str(iqr(lin.ecc(chan_sel)-DoG.ecc(chan_sel))) '\n'])

if SaveFigs
    saveas(f_neg3,fullfile(figfld,['EPHYS_NEG-PRF3_' fb{fidx} '.png']));
end
if CloseFigs; close(f_neg3); end

%% Value of exponential parameter for CSS across MUA ======================
% Run the analysis for MRI first because we're going to want to do a
% cross-signal comparison.

RTHRES = 25;

figure;
exptvals_MUA = tMUA.expt(strcmp(tMUA.Model,'css_ephys_cv1') & ...
    strcmp(tMUA.SigType,'MUA') & tMUA.R2 > RTHRES );
fprintf(['Mean expt: ' num2str(nanmean(exptvals_MUA)) ...
    ', Std ' num2str(nanstd(exptvals_MUA)) ...
    ', Median ' num2str(median(exptvals_MUA)) '\n']);
histogram(exptvals_MUA,0:.05:2)




mm = [exptv(1).roi{1};exptv(2).roi{1}];
mm = [mm ones(size(mm))];
ee = [exptvals_MUA 2*ones(size(exptvals_MUA))];

% Wilcoxon 1-tailed < 1
[p,h,stats] = signrank(ee(:,1),1,'tail','left');
fprintf(['Wilcoxon EXPT < 1 : z = ' num2str(stats.zval) ', p = ' num2str(p) '\n']);
% Mann-Whitney test
[p,h,stats] = ranksum(mm(:,1),ee(:,1));
fprintf(['Mann-Whitney U EXPT MRI vs MUA : z = ' num2str(stats.zval) ', p = ' num2str(p) '\n']);

%% Value of exponential parameter for CSS across LFP ======================
RTHRES = 25;
sig=unique(tLFP.SigType);
lfp_order = [3 1 2 5 4];
spn=1; ll=[]; 
figure;
for fb=lfp_order
    exptvals_LFP = tLFP.expt(strcmp(tLFP.Model,'css_ephys_cv1') & ...
        strcmp(tLFP.SigType,sig{fb}) & tLFP.R2 > RTHRES );
    fprintf(['Mean expt: ' num2str(nanmean(exptvals_LFP)) ...
        ', Std ' num2str(nanstd(exptvals_LFP)) ...
        ', Median ' num2str(median(exptvals_LFP)) ...
        ', IQR :' num2str(iqr(exptvals_LFP)) '\n']);
    
    subplot(1,5,spn);histogram(exptvals_LFP,0:.05:2)
    spn=spn+1;
%     mm = [exptv(1).roi{xval};exptv(2).roi{xval}];
%     mm = [mm ones(size(mm))]; % >>> MRI V1
%     ee = [exptvals_MUA 2*ones(size(exptvals_MUA))]; % >>> MUA
    ll = [ll; exptvals_LFP (2+spn)*ones(size(exptvals_LFP))];
    % Wilcoxon 1-tailed < 1
    [p,h,stats] = signrank(exptvals_LFP,1,'tail','left');
    fprintf([sig{fb} ', Wilcoxon EXPT < 1 : z = ' num2str(stats.zval) ', p = ' num2str(p) '\n']);
    % Mann-Whitney test
    [p,h,stats] = ranksum(mm(:,1),exptvals_LFP);
    fprintf(['Mann-Whitney U EXPT MRI vs LFP : z = ' num2str(stats.zval) ', p = ' num2str(p) '\n']);
end
[p,tbl,stats] = kruskalwallis(ll(:,1), ll(:,2));
[c,m,h,gnames] = multcompare(stats);

%% Manuscript comparison exponential parameter
RTHRES = 25;
clear exptvals_MUA exptvals_LFP exptvals_MRI kw
rois2 = [1 4];
f=figure;
for r=1:length(rois2)
    kw{r}=[];
    exptvals_MRI{r} = [exptv(1).roi{rois2(r)};exptv(2).roi{rois2(r)}];  
    kw{r}=[kw{r}; exptvals_MRI{r} ones(size(exptvals_MRI{r}))];
    exptvals_MUA{r} = tMUA.expt(strcmp(tMUA.Model,'css_ephys_cv1') & ...
        strcmp(tMUA.SigType,'MUA') & tMUA.Area==rois2(r) & tMUA.R2 > RTHRES );
    kw{r}=[kw{r}; exptvals_MUA{r} 2*ones(size(exptvals_MUA{r}))];
    sig=unique(tLFP.SigType);
    lfp_order = [3 1 2 5 4]; spn=1;
    
    cl=2;
    for fb=lfp_order
        cl=cl+1;
        exptvals_LFP{r,spn} = tLFP.expt(strcmp(tLFP.Model,'css_ephys_cv1') & ...
            strcmp(tLFP.SigType,sig{fb}) & tLFP.Area==rois2(r) & tLFP.R2 > RTHRES );
        if size(exptvals_LFP{r,spn},1)>3
            kw{r}=[kw{r}; exptvals_LFP{r,spn} cl*ones(size(exptvals_LFP{r,spn}))];
        end
        spn=spn+1;
    end
    
    figure(f);
    subplot(1,2,r); hold on;
    bar(1:7,[mean(exptvals_MRI{r}) ...
        mean(exptvals_MUA{r}) ...
        mean(exptvals_LFP{r,1}) ...
        mean(exptvals_LFP{r,2}) ...
        mean(exptvals_LFP{r,3}) ...
        mean(exptvals_LFP{r,4}) ...
        mean(exptvals_LFP{r,5})]);
    errorbar(1:7,[mean(exptvals_MRI{r}) ...
        mean(exptvals_MUA{r}) ...
        mean(exptvals_LFP{r,1}) ...
        mean(exptvals_LFP{r,2}) ...
        mean(exptvals_LFP{r,3}) ...
        mean(exptvals_LFP{r,4}) ...
        mean(exptvals_LFP{r,5})],...
        [std(exptvals_MRI{r}) ...
        std(exptvals_MUA{r}) ...
        std(exptvals_LFP{r,1}) ...
        std(exptvals_LFP{r,2}) ...
        std(exptvals_LFP{r,3}) ...
        std(exptvals_LFP{r,4}) ...
        std(exptvals_LFP{r,5})],'ko','Linestyle','none')
    set(gca,'xlim',[.5 7.5],'ylim',[0 1.1])
    
    
    
    fprintf('-- Compare < 1 --\n');
    fprintf(['AREA V' num2str(rois2(r)) '\n']);
    
    [p,h,stats] = signrank(exptvals_MRI{r},1,'tail','left');
    fprintf(['MRI, Wilcoxon EXPT < 1 : z = ' num2str(stats.zval) ', p = ' num2str(p) '\n']);
    [p,h,stats] = signrank(exptvals_MUA{r},1,'tail','left');
    fprintf(['MUA, Wilcoxon EXPT < 1 : z = ' num2str(stats.zval) ', p = ' num2str(p) '\n']);
    for i=1:5
        if ~isempty(exptvals_LFP{r,i})
            [p,h,stats] = signrank(exptvals_LFP{r,i},1,'tail','left');
            if isfield(stats,'zval')
                fprintf(['LFP-' num2str(i) ', Wilcoxon EXPT < 1 : z = ' num2str(stats.zval) ', p = ' num2str(p) '\n']);
            end
        end
    end
    
    fprintf('-- Compare across signals --\n');
    fprintf(['AREA V' num2str(rois2(r)) '\n']);
    [expcss(r).p,expcss(r).tbl,expcss(r).stats] = kruskalwallis(kw{r}(:,1), kw{r}(:,2));
    [expcss(r).c,expcss(r).m,expcss(r).h,expcss(r).gnames] = multcompare(expcss(r).stats);
end
close(f)





%% Manuscript comparison of location and size across ephys channels -------
RTH=25; SNRTH = 3;

% forget about theta here, it doesn't have good pRF fits

% CSS -----
modidx = 3; % 2 = U-LIN, 3 = CSS
fprintf(['MODEL ' ephys_MOD{modidx} '\n']);
% 1 = MUA, 2 = ClasRF, 3 = Theta, 4 = Alpha, 5 = Beta, 6 = Gamma-low, 7 =
% Gamma-high

SigCompNames = {};
SigComp_nElec =[];
SigCompDist = [];
SigCompDist_columns = {'mean','std','median','iqr'};
SigCompSz = [];
SigCompSz_columns = {'mean_nS1','std_nS1','median_nS1','iqr_nS1',...
    'mean_nS2','std_nS2','median_nS2','iqr_nS2',...
    'mean_nS2/nS1','std_nS2/nS1','median_nS2/nS1','iqr_nS2/nS1'};

rown=1;
for sigidx1 = 1:7
    for sigidx2 = 1:7 
        SigCompNames{rown,1} = PRF_EST(modidx,sigidx1).sig;
        SigCompNames{rown,2} = PRF_EST(modidx,sigidx2).sig;
        
        % threshold
        if sigidx1 == 2
            elec = PRF_EST(modidx,sigidx1).R2 > SNRTH & PRF_EST(modidx,sigidx2).R2 > RTH;
        elseif sigidx2 == 2
            elec = PRF_EST(modidx,sigidx1).R2 > RTH & PRF_EST(modidx,sigidx2).R2 > SNRTH;
        else
            elec = PRF_EST(modidx,sigidx1).R2 > RTH & PRF_EST(modidx,sigidx2).R2 > RTH;
        end

        % calculate distance
        DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
            (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
        SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
        
        % calculate normalized size
        nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
        nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
        SigCompSz = [SigCompSz;...
            nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
            nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
            nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
        
        SigComp_nElec = [SigComp_nElec;sum(elec)];
        
        rown = rown+1;
    end
end

[sorted_names,sidx] = sortrows(SigCompNames);
sorted_n = SigComp_nElec(sidx,:);
sorted_dist = SigCompDist(sidx,:);
sorted_sz = SigCompSz(sidx,:);

uSig = unique(SigCompNames);
nSig = size(unique(SigCompNames),1);

DistMat_median = zeros(nSig);
DistMat_iqr = zeros(nSig);
DistMat_n = zeros(nSig);
SzMat_median = zeros(nSig);
SzMat_iqr = zeros(nSig);

for i=1:nSig
    ii=((i-1)*nSig)+1;
    DistMat_median(:,i) = sorted_dist(ii:ii+nSig-1,3);
    DistMat_iqr(:,i) = sorted_dist(ii:ii+nSig-1,4)./2;
    DistMat_n(:,i) = sorted_n(ii:ii+nSig-1);
    SzMat_median(:,i) = sorted_sz(ii:ii+nSig-1,11);
    SzMat_iqr(:,i) = sorted_sz(ii:ii+nSig-1,12)./2;
end


% re-order for plot
order=[4 3 5 1 2 7 6];
DistMat_median = DistMat_median(order,:); 
DistMat_median = DistMat_median(:,order); 
DistMat_iqr = DistMat_iqr(order,:); 
DistMat_iqr = DistMat_iqr(:,order); 
SzMat_median = SzMat_median(order,:); 
SzMat_median = SzMat_median(:,order); 
SzMat_iqr = SzMat_iqr(order,:); 
SzMat_iqr = SzMat_iqr(:,order); 
DistMat_n = DistMat_n(order,:);
DistMat_n = DistMat_n(:,order); 
uSig = uSig(order);

%
fcss=figure;
set(fcss,'Position',[10 10 2200 1200],'Renderer','painters');
%colormap(viridis)
colormap(brewermap([],'RdBu'));

subplot(2,3,1);
imagesc(DistMat_median)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 1.5]);colorbar; 
title('pRF distance median');

subplot(2,3,2);
imagesc(DistMat_iqr)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 0.8]);colorbar; 
title('pRF distance iqr');

subplot(2,3,3);
imagesc(DistMat_n)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45,'ColorScale','log');
colorbar; 
caxis([1 1700]) 
title('n');

subplot(2,3,4);
imagesc(SzMat_median)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 2]);colorbar; 
title('pRF relative sz median');

subplot(2,3,5);
imagesc(SzMat_iqr)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 2]);colorbar; 
title('pRF relative sz iqr');

subplot(2,3,6);
nMask = DistMat_n >= 10;
imagesc(nMask)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45);
colorbar; 
caxis([0 1]) 
title('mask n > 10');

sgtitle(['MODEL ' ephys_MOD{modidx}], 'interpreter','none')



figure;
errorbar(1:7,SzMat_median(:,2),SzMat_iqr(:,2));

%%
% U-LIN -----
modidx = 2; % 2 = U-LIN, 3 = CSS
fprintf(['MODEL ' ephys_MOD{modidx} '\n']);
% 1 = MUA, 2 = ClasRF, 3 = Theta, 4 = Alpha, 5 = Beta, 6 = Gamma-low, 7 =
% Gamma-high

SigCompNames = {};
SigComp_nElec =[];
SigCompDist = [];
SigCompDist_columns = {'mean','std','median','iqr'};
SigCompSz = [];
SigCompSz_columns = {'mean_nS1','std_nS1','median_nS1','iqr_nS1',...
    'mean_nS2','std_nS2','median_nS2','iqr_nS2',...
    'mean_nS2/nS1','std_nS2/nS1','median_nS2/nS1','iqr_nS2/nS1'};
rown=1;
for sigidx1 = 1:7
    for sigidx2 = 1:7              
        if (sigidx1 ==  4 || sigidx1 == 5) % split signal 1
            for lfs1 = 1:3
                if lfs1 == 1 % all
                    flag1='';
                    if (sigidx2 ==  4 || sigidx2 == 5) % split signal 2
                        for lfs2 = 1:3
                            if lfs2 == 1 % all
                                flag2 = '';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH;
                            elseif lfs2 == 2 % gain > 0
                                flag2 = '_pg';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).G > 0;
                            elseif lfs2 == 3 % gain < 0
                                flag2 = '_ng';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).G < 0;
                            end
                            DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                                (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                            SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                            nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                            nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                            SigCompSz = [SigCompSz;...
                                nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                                nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                                nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                            SigComp_nElec = [SigComp_nElec;sum(elec)];
                            
                            SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                            SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
                            
                            rown = rown+1;
                        end
                    else
                        flag2 = '';
                        elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                            PRF_EST(modidx,sigidx2).R2 > RTH;
                        DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                            (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                        SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                        nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                        nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                        SigCompSz = [SigCompSz;...
                            nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                            nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                            nanmean(nS2./nS1) std(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                        SigComp_nElec = [SigComp_nElec;sum(elec)];
                        
                        SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                        SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
 
                        rown = rown+1;
                        
                    end
                elseif lfs1 == 2 % gain > 0
                    flag1 = '_pg';
                    if (sigidx2 ==  4 || sigidx2 == 5) % split signal 2
                        for lfs2 = 1:3
                            if lfs2 == 1
                                flag2='';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx1).G > 0 & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH;
                            elseif lfs2 == 2 % gain > 0
                                flag2='_pg';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx1).G > 0 & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).G > 0;
                            elseif lfs2 == 3 % gain < 0
                                flag2='_ng';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx1).G > 0 & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).G < 0;
                            end
                            DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                                (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                            SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                            nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                            nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                            SigCompSz = [SigCompSz;...
                                nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                                nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                                nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                            SigComp_nElec = [SigComp_nElec;sum(elec)];
                            
                            SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                            SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
                        
                            rown = rown+1;
                        end
                    else
                        flag2='';
                        elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                            PRF_EST(modidx,sigidx1).G > 0 & ...
                            PRF_EST(modidx,sigidx2).R2 > RTH;
                        DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                            (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                        SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                        nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                        nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                        SigCompSz = [SigCompSz;...
                            nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                            nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                            nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                        SigComp_nElec = [SigComp_nElec;sum(elec)];
                        
                        SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                        SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
                        
                        rown = rown+1;
                    end
                elseif lfs1 == 3 % gain < 0
                    flag1='_ng';
                    if (sigidx2 ==  4 || sigidx2 == 5) % split signal 2
                        for lfs2 = 1:3
                            if lfs2 == 1
                                flag2='';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx1).G < 0 & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH;
                            elseif lfs2 == 2 % gain > 0
                                flag2='_pg';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx1).G < 0 & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).G > 0;
                            elseif lfs2 == 3 % gain < 0
                                flag2='_ng';
                                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx1).G < 0 & ...
                                    PRF_EST(modidx,sigidx2).R2 > RTH & ...
                                    PRF_EST(modidx,sigidx2).G < 0;
                            end
                            DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                                (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                            SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                            nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                            nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                            SigCompSz = [SigCompSz;...
                                nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                                nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                                nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                            SigComp_nElec = [SigComp_nElec;sum(elec)];
                            
                            SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                            SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
                        
                            rown = rown+1;
                        end
                    else
                        flag2='';
                        elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                            PRF_EST(modidx,sigidx1).G < 0 & ...
                            PRF_EST(modidx,sigidx2).R2 > RTH;
                        DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                            (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                        SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                        nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                        nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                        SigCompSz = [SigCompSz;...
                            nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                            nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                            nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                        SigComp_nElec = [SigComp_nElec;sum(elec)];
                        
                        SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                        SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
                        
                        rown = rown+1;
                    end
                end
            end
        else %do not split signal 1
            flag1='';
            if (sigidx2 ==  4 || sigidx2 == 5) % split signal 2
                for lfs2 = 1:3
                    if lfs2 == 1
                        flag2='';
                        elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                            PRF_EST(modidx,sigidx2).R2 > RTH;
                    elseif lfs2 == 2 % gain > 0
                        flag2='_pg';
                        elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                            PRF_EST(modidx,sigidx2).R2 >  RTH & ...
                            PRF_EST(modidx,sigidx2).G > 0;
                    elseif lfs2 == 3 % gain < 0
                        flag2='_ng';
                        elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                            PRF_EST(modidx,sigidx2).R2 > RTH & ...
                            PRF_EST(modidx,sigidx2).G < 0;
                    end
                    DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                        (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                    SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                    nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                    nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                    SigCompSz = [SigCompSz;...
                        nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                        nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                        nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                    SigComp_nElec = [SigComp_nElec;sum(elec)];
                    
                    SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                    SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
                        
                    rown = rown+1;
                end
            else
                flag2='';
                elec = PRF_EST(modidx,sigidx1).R2 > RTH & ...
                    PRF_EST(modidx,sigidx2).R2 > RTH;
                DIST = sqrt((PRF_EST(modidx,sigidx1).X(elec) - PRF_EST(modidx,sigidx2).X(elec)).^2 + ...
                    (PRF_EST(modidx,sigidx1).Y(elec) - PRF_EST(modidx,sigidx2).Y(elec)).^2);
                SigCompDist = [SigCompDist; mean(DIST) std(DIST) median(DIST) iqr(DIST)];
                nS1 = PRF_EST(modidx,sigidx1).S(elec)./PRF_EST(modidx,1).S(elec);
                nS2 = PRF_EST(modidx,sigidx2).S(elec)./PRF_EST(modidx,1).S(elec);
                SigCompSz = [SigCompSz;...
                    nanmean(nS1) nanstd(nS1) nanmedian(nS1) iqr(nS1) ...
                    nanmean(nS2) nanstd(nS2) nanmedian(nS2) iqr(nS2) ...
                    nanmean(nS2./nS1) nanstd(nS2./nS1) nanmedian(nS2./nS1) iqr(nS2./nS1)];
                SigComp_nElec = [SigComp_nElec;sum(elec)];
                
                SigCompNames{rown,1} = [PRF_EST(modidx,sigidx1).sig flag1]; 
                SigCompNames{rown,2} = [PRF_EST(modidx,sigidx2).sig flag2];
                        
                rown = rown+1;
            end
        end      
    end
end

[sorted_names,sidx] = sortrows(SigCompNames);
sorted_n = SigComp_nElec(sidx,:);
sorted_dist = SigCompDist(sidx,:);
sorted_sz = SigCompSz(sidx,:);

uSig = unique(SigCompNames);
nSig = size(unique(SigCompNames),1);

DistMat_median = zeros(nSig);
DistMat_iqr = zeros(nSig);
DistMat_n = zeros(nSig);
SzMat_median = zeros(nSig);
SzMat_iqr = zeros(nSig);

for i=1:nSig
    ii=((i-1)*nSig)+1;
    DistMat_median(:,i) = sorted_dist(ii:ii+nSig-1,3);
    DistMat_iqr(:,i) = sorted_dist(ii:ii+nSig-1,4)./2;
    DistMat_n(:,i) = sorted_n(ii:ii+nSig-1);
    SzMat_median(:,i) = sorted_sz(ii:ii+nSig-1,11);
    SzMat_iqr(:,i) = sorted_sz(ii:ii+nSig-1,12)./2;
end


% re-order for plot
order=[8 7 9 1:6 11 10];
DistMat_median = DistMat_median(order,:); 
DistMat_median = DistMat_median(:,order); 
DistMat_iqr = DistMat_iqr(order,:); 
DistMat_iqr = DistMat_iqr(:,order); 
SzMat_median = SzMat_median(order,:); 
SzMat_median = SzMat_median(:,order); 
SzMat_iqr = SzMat_iqr(order,:); 
SzMat_iqr = SzMat_iqr(:,order); 
DistMat_n = DistMat_n(order,:);
DistMat_n = DistMat_n(:,order); 
uSig = uSig(order);



%
fulin=figure;
set(fulin,'Position',[10 10 2200 1200],'Renderer','painters');
colormap(brewermap([],'RdBu'));
%colormap(viridis)

subplot(2,3,1);
imagesc(DistMat_median)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 3]);colorbar; 
title('pRF distance median');

subplot(2,3,2);
imagesc(DistMat_iqr)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 3]);colorbar; 
title('pRF distance iqr');

subplot(2,3,3);
imagesc(DistMat_n)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45,'ColorScale','log');
colorbar; 
caxis([1 1700]) 
title('n');

subplot(2,3,4);
imagesc(SzMat_median)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 2]);colorbar; 
title('pRF relative sz median');

subplot(2,3,5);
imagesc(SzMat_iqr)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45)
caxis([0 2]);colorbar; 
title('pRF relative sz iqr');

subplot(2,3,6);
nMask = DistMat_n >= 10;
imagesc(nMask)
set(gca,'TickDir','out','xtick',1:11,'xticklabels',uSig, 'yticklabels',uSig,...
    'XTickLabelRotation',45);
colorbar; 
caxis([0 1]) 
title('mask n > 10');

sgtitle(['MODEL ' ephys_MOD{modidx}], 'interpreter','none')

figure;
errorbar(1:11,SzMat_median(:,2),SzMat_iqr(:,2));