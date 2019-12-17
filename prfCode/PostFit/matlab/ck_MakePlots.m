% ck_MakePlots
% Takes the pre-processed fitting result tables and creates comparison plots

%% Paths ==================================================================
BaseFld = pwd;
DS='ORG';
ResFld = ...
    ['/Users/chris/Documents/MRI_ANALYSIS/NHP-PRF/FitResults/MultiModal/'...
    DS '/cv1'];
T='Tables_max';

% add colorbrewer
addpath('/Users/chris/Dropbox/MATLAB_NONGIT/TOOLBOX/BrewerMap')
def_cmap = 'Spectral';

%% Load ===================================================================
fprintf('Loading results table. Please be patient, this will take a while..\n');
tic; load(fullfile(ResFld,T)); t_load=toc;
fprintf(['Loading took ' num2str(t_load) ' s\n']);

%% Split by model =========================================================
clear T
m = unique(tMRI_max.Model); 
for mi = 1:length(m)
    M = tMRI_max(strcmp(tMRI_max.Model,m{mi}),:);
    T(mi).mod = M;
    T(mi).name = m{mi};
    modidx.(m{mi}) = mi;
end

%% scatter plots & differences R2 =========================================
roi={'V1','V2_merged','V3_merged','V4_merged','MT','MST','TEO','LIP_merged'};
roilabels={'V1','V2','V3','V4','MT','MST','TEO','LIP'};
MRI_MODELS={...
    'linhrf_cv1_mhrf','linhrf_cv1_dhrf';...
    'csshrf_cv1_mhrf','csshrf_cv1_dhrf';...
    'doghrf_cv1_mhrf','doghrf_cv1_dhrf';...
    };

f=figure;
set(f,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
set(f,'Position',[100 100 1600 1000]);
s_R2 = T(modidx.linhrf_cv1_mhrf).mod.R2>0; % allows selection but keep at 0

subplot(2,3,1); hold on;
plot([0 100],[0 100],'k','Linewidth',2);
for r=1:length(roi)
    scatter(T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),100,'Marker','.');
end
title('linear vs css'); xlabel('linear');ylabel('css');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

subplot(2,3,2); hold on;
plot([0 100],[0 100],'k','Linewidth',2); 
for r=1:length(roi)
    scatter(T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r})),100,'Marker','.');
end
title('linear vs dog'); xlabel('linear');ylabel('dog');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

subplot(2,3,3); hold on;
plot([0 100],[0 100],'k','Linewidth',2); 
for r=1:length(roi)
    scatter(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r})),100,'Marker','.');
end
title('css vs dog'); xlabel('css');ylabel('dog');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
legend(['XY' roilabels],'interpreter','none','Location','SouthEast');

% diff distrutions plots -----
diffmat{1}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(modidx.csshrf_cv1_mhrf).mod.(roi{r}));
    
    m = mean(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})));
    sd = std(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(modidx.csshrf_cv1_mhrf).mod.(roi{r})));
    diffmat{1} = [diffmat{1}; m sd se];
end

diffmat{2}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(modidx.doghrf_cv1_mhrf).mod.(roi{r}));
    
    m = mean(T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})));
    sd = std(T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(modidx.doghrf_cv1_mhrf).mod.(roi{r})));
    diffmat{2} = [diffmat{2}; m sd se];   
end

diffmat{3}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(modidx.csshrf_cv1_mhrf).mod.(roi{r}));
    
    m = mean(T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})));
    sd = std(T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(modidx.csshrf_cv1_mhrf).mod.(roi{r})));
    diffmat{3} = [diffmat{3}; m sd se];
end

TITLES={...
    'CSS - LINEAR',...
    'DOG - LINEAR',...
    'DOG - CSS'};
for cc = 1:length(diffmat)
    subplot(2,3,3+cc); hold on;

    for xval=1:length(diffmat{cc})
        bar(xval,diffmat{cc}(xval,1));
    end
    for xval=1:length(diffmat{cc})
        errorbar(xval,diffmat{cc}(xval,1),diffmat{cc}(xval,3),...
            'k-','Linestyle','none');
    end
    set(gca,'xticklabels',[],'ylim',[-1.5 2]);
    xlabel('ROI'); ylabel('Diff R2');
    title(TITLES{cc});
    if cc==3
        legend(roilabels,'interpreter','none','Location','NorthEast');
    end
end

%% scatter plot HRF & differences =========================================
f2=figure;
set(f2,'Position',[100 100 1000 1200]);
set(f2,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
s_R2 = T(modidx.linhrf_cv1_mhrf).mod.R2>0;


for mm = 1:size(MRI_MODELS,1)

    subplot(3,2,(mm-1)*2 +1); hold on;
    plot([0 100],[0 100],'k','LineWidth',2);
    for r=1:length(roi)
        scatter(T(modidx.(MRI_MODELS{mm,1})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r})),...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,2})).mod.(roi{r})),100,'Marker','.');
    end
    title(['mHRF vs cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
    xlabel('Monkey HRF R2'); ylabel('Canonical HRF R2');
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

    diffmat2{1}=[];
    for r=1:length(roi)
        [n,x] = hist(...
            T(modidx.(MRI_MODELS{mm,1})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r}))-...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,2})).mod.(roi{r})),100);
        f = n./sum(s_R2 & T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r}));
    
        m = mean(T(modidx.(MRI_MODELS{mm,1})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r}))-...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,2})).mod.(roi{r})));
        sd = std(T(modidx.(MRI_MODELS{mm,1})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r}))-...
            T(modidx.(MRI_MODELS{mm,2})).mod.R2(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,2})).mod.(roi{r})));
        se = sd ./ sqrt(sum(s_R2 & T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r})));
        diffmat2{1} = [diffmat2{1}; m sd se];
    end

    subplot(3,2,(mm-1)*2 +2); hold on
    for xval=1:length(diffmat2{1})
        bar(xval,diffmat2{1}(xval,1));
    end
    for xval=1:length(diffmat2{1})
        errorbar(xval,diffmat2{1}(xval,1),diffmat2{1}(xval,3),...
        'k-','Linestyle','none')
    end
    set(gca,'xticklabels',[],'ylim',[-0.65 1]);
    xlabel('ROI'); ylabel('Diff R2');
    title(['mHRF - cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
    legend(roilabels,'interpreter','none','Location','NorthEast');
end


%% rf size depending on HRF ===============================================
f3=figure;
set(f3,'Position',[100 100 1000 1200]);
set(f3,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));

R2th=10;
s_R2 = T(modidx.csshrf_cv1_mhrf).mod.R2 > R2th;

for mm = 1:size(MRI_MODELS,1)
    
    sp=subplot(3,2,(mm-1)*2 +1);hold on;
    plot([0 100],[0 100],'k','LineWidth',2);
    for r=1:length(roi)
        scatter(T(modidx.(MRI_MODELS{mm,1})).mod.rfs(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r})),...
            T(modidx.(MRI_MODELS{mm,2})).mod.rfs(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,2})).mod.(roi{r})),100,'Marker','.');
    end
    title(['mHRF vs cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
    xlabel('Monkey HRF sz');ylabel('Canonical HRF sz');
    set(gca, 'Box','off', 'xlim', [0 15], 'ylim',[0 15]);
    text('Parent',sp,'Position',[1 13], ...
        'String',['R2th: ' num2str(R2th)],...
        'Fontsize',12, 'Fontweight','bold')

    
    diffmat2{1}=[];
    for r=1:length(roi)
        [n,x] = hist(...
            T(modidx.(MRI_MODELS{mm,1})).mod.rfs(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r}))-...
            T(modidx.(MRI_MODELS{mm,2})).mod.rfs(s_R2 & ...
            T(modidx.(MRI_MODELS{mm,2})).mod.(roi{r})),...
            100);
        f = n./sum(s_R2 & T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r}));
        
        dsz = T(modidx.(MRI_MODELS{mm,1})).mod.rfs-T(modidx.(MRI_MODELS{mm,2})).mod.rfs;
        dsz = dsz(s_R2 & T(modidx.(MRI_MODELS{mm,1})).mod.(roi{r}));
        dsz = dsz(isfinite(dsz));
        
        m = mean(dsz);
        sd = std(dsz);
        se = sd ./ sqrt(length(dsz));
        diffmat2{1} = [diffmat2{1}; m sd se];
    end
    
    
    subplot(3,2,(mm-1)*2 +2); hold on
    for xval=1:length(diffmat2{1})
        bar(xval,diffmat2{1}(xval,1));
    end
    for xval=1:length(diffmat2{1})
        errorbar(xval,diffmat2{1}(xval,1),diffmat2{1}(xval,3),...
            'k-','Linestyle','none')
    end
    set(gca,'xticklabels',[],'ylim',[-0.5 1.5]);
    xlabel('ROI'); ylabel('Diff pRF size');
    title(['mHRF - cHRF (' MRI_MODELS{mm,1} ')'],'interpreter','none'); 
    legend(roilabels,'interpreter','none','Location','NorthEast');
end
%% ECC vs PRF Size ========================================================
f4=figure;
set(f4,'Position',[100 100 1000 1200]);
set(f4,'DefaultAxesColorOrder',brewermap(length(roi),def_cmap));
Rth=15; 
MRI_MODEL={'linhrf_cv1_mhrf','csshrf_cv1_mhrf','doghrf_cv1_mhrf'};

for m=1:length(MRI_MODEL)
    s_R2 = T(modidx.(MRI_MODEL{m})).mod.R2 > Rth;
    EccBin = 0.5:1:30.5;
    
    subplot(3,2,(m-1)*2 +1);hold on;
    for r=1:length(roi)
        ES{r}=[];
        scatter(T(modidx.(MRI_MODEL{m})).mod.ecc(s_R2 & ...
            T(modidx.(MRI_MODEL{m})).mod.(roi{r})),...
            T(modidx.(MRI_MODEL{m})).mod.rfs(s_R2 & ...
            T(modidx.(MRI_MODEL{m})).mod.(roi{r})),100,'Marker','.');
        for b=1:length(EccBin)
            bb=[EccBin(b)-0.5 EccBin(b)+0.5];
            PSZ=T(modidx.(MRI_MODEL{m})).mod.rfs(s_R2 & ...
                T(modidx.(MRI_MODEL{m})).mod.(roi{r}));
            ECC=T(modidx.(MRI_MODEL{m})).mod.ecc(s_R2 & ...
                T(modidx.(MRI_MODEL{m})).mod.(roi{r}));
            ES{r}=[ES{r}; EccBin(b) median(PSZ(ECC>=bb(1) & ECC<=bb(2)))];
        end
    end
    title(['Ecc vs pRF size [' MRI_MODEL{m} ', R>' num2str(Rth) ']'],...
        'interpreter','none');
    xlabel('Eccentricity');ylabel('pRF size');
    set(gca, 'Box','off', 'xlim', [0 10], 'ylim',[0 10]);
    
    subplot(3,2,(m-1)*2 +2);hold on;
    for r=1:length(roi)
        h=plot(ES{r}(:,1),ES{r}(:,2),'o');
        set(h,'MarkerSize',6,'markerfacecolor', get(h, 'color'));
    end
    title(['Ecc vs pRF size [' MRI_MODEL{m} ', R>' num2str(Rth) ']'],...
        'interpreter','none'); xlabel('Eccentricity');ylabel('pRF size');
    set(gca, 'Box','off', 'xlim', [0 10], 'ylim',[0 10]);
    legend(roilabels,'interpreter','none','Location','NorthWest');
end

%% Ephys VFC ==============================================================
% MUA
s=tMUA_max.R2>25;
szm=tMUA_max.rfs<200;

model='css_ephys_cv1';

fm=figure;
set(fm,'Position',[100 100 800 800]);
%[cmap,~,~] = brewermap(length(unique(tMUA_max.Array)),'Dark2');
%set(fm,'DefaultAxesColorOrder',cmap);

subplot(2,2,1);hold on;

szm=tMUA_max.rfs<5;
m=strcmp(tMUA_max.Monkey,'lick') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==1;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v),...
        tMUA_max.Y(s & m & a & v),'Marker','*' )
end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Lick V1 MUA')

subplot(2,2,3);hold on;
szm=tMUA_max.rfs<5;
m=strcmp(tMUA_max.Monkey,'lick') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==4;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v),...
        tMUA_max.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Lick V4 MUA')

subplot(2,2,2);hold on;
m=strcmp(tMUA_max.Monkey,'aston') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==1;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v),...
        tMUA_max.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Aston V1 MUA')

subplot(2,2,4);hold on;
m=strcmp(tMUA_max.Monkey,'aston') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==4;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v),...
        tMUA_max.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Aston V4 MUA')

%% Ephys VFC ==============================================================
% LFP Low Gamma
s=tLFP_max.R2>25;

model='css_ephys_cv1';
b = 'lGamma';
fm=figure;
set(fm,'Position',[100 100 800 800]);

subplot(2,2,1);hold on;
szm=tLFP_max.rfs<50;
m=strcmp(tLFP_max.Monkey,'lick') & ...
    strcmp(tLFP_max.Model, model) & ...
    strcmp(tLFP_max.SigType, b);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP_max.Area==1;
for r=unique(tLFP_max.Array)'
    a=tLFP_max.Array==r;
    scatter(tLFP_max.X(s & m & a & v),...
        tLFP_max.Y(s & m & a & v),'Marker','*' )
end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Lick V1 LFP (lGAM)')

subplot(2,2,3);hold on;
szm=tLFP_max.rfs<50;
m=strcmp(tLFP_max.Monkey,'lick') & strcmp(tLFP_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP_max.Area==4;
for r=unique(tLFP_max.Array)'
    a=tLFP_max.Array==r;
    scatter(tLFP_max.X(s & m & a & v),...
        tLFP_max.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Lick V4 LFP (lGAM)')

subplot(2,2,2);hold on;
m=strcmp(tLFP_max.Monkey,'aston') & strcmp(tLFP_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP_max.Area==1;
for r=unique(tLFP_max.Array)'
    a=tLFP_max.Array==r;
    scatter(tLFP_max.X(s & m & a & v),...
        tLFP_max.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Aston V1 LFP (lGAM)')

subplot(2,2,4);hold on;
m=strcmp(tLFP_max.Monkey,'aston') & strcmp(tLFP_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP_max.Area==4;
for r=unique(tLFP_max.Array)'
    a=tLFP_max.Array==r;
    scatter(tLFP_max.X(s & m & a & v),...
        tLFP_max.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-5 10], 'ylim',[-10 5]);
title('Aston V4 LFP (lGAM)')

%% Ephys location difference & size difference  ---
rth=25;

ephys_MOD={'linear_ephys_cv1','css_ephys_cv1','dog_ephys_cv1'};
clear C R2m SZ
for m=1:length(ephys_MOD)
 
    C{m}=[];R2m{m}=[];SZ{m}=[];
    
    model=ephys_MOD{m};
    s = strcmp(tMUA_max.Model,model);
    C{m}=[C{m} tMUA_max.R2(s) tMUA_max.X(s) tMUA_max.Y(s) tMUA_max.rfs(s)];
    R2m{m}=[R2m{m} tMUA_max.R2(s)];
    SZ{m}=[SZ{m} tMUA_max.R2(s) tMUA_max.rfs(s)];
    
    s = strcmp(tMUA_max.Model,'classicRF');
    C{m}=[C{m} tMUA_max.X(s)./668.745 tMUA_max.Y(s)./668.745 tMUA_max.rfs(s)./2];
    
    s = strcmp(tLFP_max.Model,model);
    sig=unique(tLFP_max.SigType);
    lfp_order = [3 1 2 5 4];
    for i=lfp_order
        b=strcmp(tLFP_max.SigType,sig{i});
        C{m}=[C{m} tLFP_max.R2(s & b) tLFP_max.X(s & b) tLFP_max.Y(s & b) tLFP_max.rfs(s & b)];
        R2m{m}=[R2m{m} tLFP_max.R2(s & b)];
        SZ{m}=[SZ{m} tLFP_max.R2(s & b) tLFP_max.rfs(s & b)];
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
%% MUA model comparison ===================================================
m=unique(tMUA_max.Model);
R2=[];
for i=1:length(m)
    R2 = [R2 tMUA_max.R2(strcmp(tMUA_max.Model,m{i}))];
end

v1=tLFP_max.Area(strcmp(tMUA_max.Model,m{1}))==1;
v4=tLFP_max.Area(strcmp(tMUA_max.Model,m{1}))==4;

f=figure;
set(f,'Position',[100 100 1200 1200]);
for row=1:3
    for column=1:3
        subplot(3,3,((row-1)*3)+column); hold on;
        plot([0 100],[0 100],'k');
        scatter(R2(v1,row+1), R2(v1,column+1),60,'Marker','.',...
            'MarkerEdgeColor',[.3 .3 .3]);
        scatter(R2(v4,row+1), R2(v4,column+1),60,'Marker','.',...
            'MarkerEdgeColor',[.3 .8 .3]);
        set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
        xlabel(m{row+1},'interpreter','none'); 
        ylabel(m{column+1},'interpreter','none');
        title('MUA');
        if row==1 && column==1
            legend({'','V1','V4'},'location','NorthWest');
        end
    end
end

%% LFP model comparison ===================================================
f=figure; 
set(f,'Position',[100 100 1200 1200]);

sig=unique(tLFP_max.SigType);
lfp_order = [3 1 2 5 4];
spn=1;
for fb=lfp_order
    m=unique(tLFP_max.Model);
    R2=[];
    for i=1:length(m)-1
        R2 = [R2 tLFP_max.R2(...
            strcmp(tLFP_max.Model,m{i}) & ...
            strcmp(tLFP_max.SigType,sig{fb}))];
    end
    
    subplot(length(sig),3,spn); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(v1,3), R2(v1,1),60,'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    scatter(R2(v4,3), R2(v4,1),60,'Marker','.',...
        'MarkerEdgeColor',[.3 .8 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    xlabel(m{3},'interpreter','none');ylabel(m{1},'interpreter','none')
    title(sig{fb})
    if spn==1
        legend({'','V1','V4'},'location','SouthEast');
    end
    spn=spn+1;
    
    subplot(length(sig),3,spn); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(v1,3), R2(v1,2),60,'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    scatter(R2(v4,3), R2(v4,2),60,'Marker','.',...
        'MarkerEdgeColor',[.3 .8 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    xlabel(m{3},'interpreter','none');ylabel(m{2},'interpreter','none')
    title(sig{fb})
    spn=spn+1;
    
    subplot(length(sig),3,spn); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(v1,1), R2(v1,2),60,'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    scatter(R2(v4,1), R2(v4,2),60,'Marker','.',...
        'MarkerEdgeColor',[.3 .8 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    xlabel(m{1},'interpreter','none');ylabel(m{2},'interpreter','none')
    title(sig{fb})
    spn=spn+1;
      
end

%% R2 for different ephys signals =========================================
RR=[];
ephys_MOD={'linear_ephys_cv1','css_ephys_cv1','dog_ephys_cv1'};

for m=1:length(ephys_MOD)
    model=ephys_MOD{m};
    s = strcmp(tMUA_max.Model,model);
    RR=[RR tMUA_max.R2(s)];
    
    s = strcmp(tLFP_max.Model,model);
    sig=unique(tLFP_max.SigType);
    lfp_order = [3 1 2 5 4];
    for i=lfp_order
        b=strcmp(tLFP_max.SigType,sig{i});
        RR=[RR tLFP_max.R2(s & b)];
    end
    LAB=['MUA';sig(lfp_order)];
    
    f=figure; set(f,'Position',[100 100 1300 1200]);
    r2th=0;
    
    c=0;d=0;
    for ref=1:6
        c=c+1;
        for fb=1:6
            d=d+1;
            s=(RR(:,ref)>r2th & RR(:,fb)>r2th);
            subplot(6,6,d); hold on; plot([0 100],[0 100],'k');
            scatter(RR(s,ref),RR(s,fb),120,[0.3 0.3 0.3],'Marker','.')
            xlabel(LAB{ref});ylabel(LAB{fb});
            title(['R2 ' model],'interpreter','none');
            set(gca,'xlim',[0 100],'ylim',[0 100]);
        end
    end

    % Distance from diagonal ==============================================
    f=figure; set(f,'Position',[100 100 1300 1200]);
    LAB=['MUA';sig(lfp_order)];

    c=0;d=0;
    for ref=1:6
        c=c+1;
        for fb=1:6
            d=d+1;
            s=(RR(:,ref)>r2th & RR(:,fb)>r2th);
            
            subplot(6,6,d); hold on;
            dRR = RR(s,fb)-RR(s,ref);
            h = histogram(dRR,-100:1:100,'FaceColor','k','FaceAlpha',1);
            YY = get(gca,'ylim');
            plot([0 0],YY,'Color',[0.5 0.5 0.5],'LineWidth',2)
            plot([mean(dRR) mean(dRR)],YY,'r','LineWidth',2)
            plot([median(dRR) median(dRR)],YY,'b','LineWidth',2)
            set(gca,'xlim',[-100 100])
            xlabel(['dRR ' model],'interpreter','none');
            ylabel('cnt','interpreter','none');
            title(['R2 ' LAB{(fb)} '-' LAB{(ref)} ],...
                'interpreter','none')
            
            if ref ==1 && fb ==1
                legend({'HIST','0','MEAN','MEDIAN'});
            end
            
        end
    end
end


%% pRF size for different ephys signals ===================================
% SZ is [ MUA_R2(1) MUA_RFS(2) 
%         THETA_R2(3) THETA_RFS(4) 
%         ALPHA_R2(5) ALPHA_RFS(6) 
%         BETA_R2(7) BETA_RFS8) 
%         LGAM_R2(9) LG_RFS(10) 
%         HGAM_R2(11) HGAM_RFS(12)]


ephys_MOD={'linear_ephys_cv1','css_ephys_cv1','dog_ephys_cv1'};

for m=1:length(ephys_MOD)
    model=ephys_MOD{m};
    LAB=['MUA';sig(lfp_order)];
    
    f=figure; set(f,'Position',[100 100 1300 1200]);
    r2th=20;
    
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
            title(['pRF size' model(1:3)],'interpreter','none')
            
            set(gca,'xlim',[0 10],'ylim',[0 10]);
            
            SS(c).nSZ = [SS(c).nSZ ; ...
                median(diff(SZ{m}(s,[ref+1 fb+1]),1,2)) ...
                std(diff(SZ{m}(s,[ref+1 fb+1]),1,2))./sqrt(sum(s)) ...
                median(  SZ{m}(s,fb+1)./SZ{m}(s,ref+1) ) ...
                std(  SZ{m}(s,fb+1)./SZ{m}(s,ref+1) )./sqrt(sum(s))];
        end
    end
    
    % Distance from diagonal ==============================================
    
    f=figure; set(f,'Position',[100 100 1300 1200]);
    r2th=20;
    
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
            xlabel(['Sz ' model(1:3)],'interpreter','none');
            ylabel('cnt','interpreter','none');
            title(['Size ' LAB{(fb+1)/2} '-' LAB{(ref+1)/2} ],...
                'interpreter','none')
            
            if ref ==1 && fb ==1
                legend({'HIST','0','MEAN','MEDIAN'});
            end
            
        end
    end
end


%% What's specific about the good DoG fits ================================
% - find these channels/voxels
% - plot their location
% - plot their size
% - plot their prf profile





%% Pattern Space-Size =====================================================
% - classify all pRFs by their location and get their size
% - correlate across modailities

