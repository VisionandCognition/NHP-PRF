% ck_MakePlots
% Takes the pre-processed fitting result tables and creates comparison plots

%% Paths ==================================================================
BaseFld = pwd;
ResFld = ...
    '/Users/chris/Documents/MRI_ANALYSIS/NHP-PRF/FitResults/MultiModal/cv1';
T='Tables_max';

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
f=figure;
set(f,'Position',[100 100 1600 1000]);
s_R2 = T(modidx.linhrf_cv1_mhrf).mod.R2>0;
roi={'V1','V2_merged','V3_merged','V4_merged','MT','MST','TEO','LIP_merged'};
roilabels={'V1','V2','V3','V4','MT','MST','TEO','LIP'};

subplot(2,3,1); hold on;
plot([0 100],[0 100],'k','Linewidth',2);
% scatter(T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2),...
%     T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2),'Marker','.',...
%     'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),'Marker','.');
end
title('linear vs css'); xlabel('linear');ylabel('css');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

subplot(2,3,2); hold on;
plot([0 100],[0 100],'k','Linewidth',2); 
% scatter(T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2),...
%     T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2),'Marker','.',...
%     'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(modidx.linhrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.linhrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r})),'Marker','.');
end
title('linear vs dog'); xlabel('linear');ylabel('dog');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

subplot(2,3,3); hold on;
plot([0 100],[0 100],'k','Linewidth',2); 
% scatter(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2),...
%     T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2),'Marker','.',...
%     'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.doghrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.doghrf_cv1_mhrf).mod.(roi{r})),'Marker','.');
end
title('css vs dog'); xlabel('css');ylabel('dog');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
legend(['XY' roi],'interpreter','none','Location','SouthEast');

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
set(f2,'Position',[100 100 1000 400]);

subplot(1,2,1); hold on;
plot([0 100],[0 100],'k','LineWidth',2);
% scatter(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2),...
%     T(modidx.csshrf_cv1_dhrf).mod.R2(s_R2),'Marker','.',...
%     'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.csshrf_cv1_dhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_dhrf).mod.(roi{r})),'Marker','.');
end
title('Monkey vs Canonical HRF'); 
xlabel('Monkey HRF R2'); ylabel('Canonical HRF R2');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

diffmat2{1}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.csshrf_cv1_dhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_dhrf).mod.(roi{r})),100);
    f = n./sum(s_R2 & T(modidx.csshrf_cv1_mhrf).mod.(roi{r}));
    
    m = mean(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.csshrf_cv1_dhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_dhrf).mod.(roi{r})));
    sd = std(T(modidx.csshrf_cv1_mhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r}))-...
        T(modidx.csshrf_cv1_dhrf).mod.R2(s_R2 & ...
        T(modidx.csshrf_cv1_dhrf).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(modidx.csshrf_cv1_mhrf).mod.(roi{r})));
    diffmat2{1} = [diffmat2{1}; m sd se];
end

subplot(1,2,2); hold on
for xval=1:length(diffmat2{1})
    bar(xval,diffmat2{1}(xval,1));
end
for xval=1:length(diffmat2{1})
    errorbar(xval,diffmat2{1}(xval,1),diffmat2{1}(xval,3),...
    'k-','Linestyle','none')
end
set(gca,'xticklabels',[],'ylim',[-0.5 1]);
xlabel('ROI'); ylabel('Diff R2');
title('Mon.HRF - Can.HRF');
legend(roilabels,'interpreter','none','Location','NorthEast');

%% rf size depending on HRF ===============================================
f3=figure;
set(f3,'Position',[100 100 1000 400]);

s_R2 = T(modidx.csshrf_cv1_mhrf).mod.R2 > 5;

MOD={'csshrf_cv1_mhrf','csshrf_cv1_dhrf'};
%MOD={'linhrf_cv1_mhrf','linhrf_cv1_dhrf'};

subplot(1,2,1);hold on;
plot([0 100],[0 100],'k','LineWidth',2);
% scatter(T(modidx.csshrf_cv1_mhrf).mod.rfs(s_R2),...
%     T(modidx.csshrf_cv1_dhrf).mod.rfs(s_R2),'Marker','.',...
%     'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(modidx.(MOD{1})).mod.rfs(s_R2 & ...
        T(modidx.(MOD{1})).mod.(roi{r})),...
        T(modidx.csshrf_cv1_dhrf).mod.rfs(s_R2 & ...
        T(modidx.csshrf_cv1_dhrf).mod.(roi{r})),'Marker','.');
end
title('Monkey HRF vs Canonical HRF'); 
xlabel('Monkey HRF sz');ylabel('Canonical HRF sz');
set(gca, 'Box','off', 'xlim', [0 15], 'ylim',[0 15]);

diffmat2{1}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(modidx.(MOD{1})).mod.rfs(s_R2 & ...
        T(modidx.(MOD{1})).mod.(roi{r}))-...
        T(modidx.(MOD{2})).mod.rfs(s_R2 & ...
        T(modidx.(MOD{2})).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(modidx.(MOD{1})).mod.(roi{r}));
    
    dsz = T(modidx.(MOD{1})).mod.rfs-T(modidx.(MOD{2})).mod.rfs;
    dsz = dsz(s_R2 & T(modidx.(MOD{1})).mod.(roi{r}));
    dsz = dsz(isfinite(dsz));
    
    m = mean(dsz);
    sd = std(dsz);
    se = sd ./ sqrt(length(dsz));
    diffmat2{1} = [diffmat2{1}; m sd se];
end


subplot(1,2,2); hold on
for xval=1:length(diffmat2{1})
    bar(xval,diffmat2{1}(xval,1));
end
for xval=1:length(diffmat2{1})
    errorbar(xval,diffmat2{1}(xval,1),diffmat2{1}(xval,3),...
    'k-','Linestyle','none')
end
set(gca,'xticklabels',[],'ylim',[-0.5 1.5]);
xlabel('ROI'); ylabel('Diff pRF size');
title('Mon.HRF - Can.HRF');
legend(roilabels,'interpreter','none','Location','NorthEast');

%% ECC vs PRF Size ========================================================
f4=figure;
subplot(1,2,1);hold on;
set(f4,'Position',[100 100 1000 400]);

s_R2 = T(modidx.csshrf_cv1_mhrf).mod.R2 > 10;
% scatter(T(modidx.csshrf_cv1_mhrf).mod.ecc(s_R2),...
%     T(modidx.csshrf_cv1_mhrf).mod.rfs(s_R2),'Marker','.',...
%     'MarkerEdgeColor',[.3 .3 .3]);
EccBin = 0.5:1:30.5;

for r=1:length(roi)
    ES{r}=[];
    scatter(T(modidx.csshrf_cv1_mhrf).mod.ecc(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),...
        T(modidx.csshrf_cv1_mhrf).mod.rfs(s_R2 & ...
        T(modidx.csshrf_cv1_mhrf).mod.(roi{r})),'Marker','.');
    for b=1:length(EccBin) 
        bb=[EccBin(b)-0.5 EccBin(b)+0.5];
        PSZ=T(modidx.csshrf_cv1_mhrf).mod.rfs(s_R2 & ...
            T(modidx.csshrf_cv1_mhrf).mod.(roi{r}));
        ECC=T(modidx.csshrf_cv1_mhrf).mod.ecc(s_R2 & ...
            T(modidx.csshrf_cv1_mhrf).mod.(roi{r}));
        ES{r}=[ES{r}; EccBin(b) mean(PSZ(ECC>=bb(1) & ECC<=bb(2)))];
    end    
end
title('Eccentricity vs pRF size'); 
xlabel('Eccentricity');ylabel('pRF size');
set(gca, 'Box','off', 'xlim', [0 15], 'ylim',[0 10]);
%legend(roilabels,'interpreter','none','Location','NorthEast');

subplot(1,2,2);hold on;
for r=1:length(roi)
    h=plot(ES{r}(:,1),ES{r}(:,2),'o');
    set(h,'MarkerSize',6,'markerfacecolor', get(h, 'color'));
end
title('Eccentricity vs pRF size'); 
xlabel('Eccentricity');ylabel('pRF size');
set(gca, 'Box','off', 'xlim', [0 15], 'ylim',[0 10]);
legend(roilabels,'interpreter','none','Location','NorthWest');

%% Ephys VFC ==============================================================
% MUA
s=tMUA_max.R2>25;
szm=tMUA_max.rfs<200;

model='css_ephys_cv1';

fm=figure;
set(fm,'Position',[100 100 800 800]);

subplot(2,2,1);hold on;
szm=tMUA_max.rfs<5;
m=strcmp(tMUA_max.Monkey,'lick') & strcmp(tMUA_max.Model,model);;
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
C=[];R2m=[];SZ=[];

model='css_ephys_cv1';
s = strcmp(tMUA_max.Model,model);
C=[C tMUA_max.R2(s) tMUA_max.X(s) tMUA_max.Y(s) tMUA_max.rfs(s)]; 
R2m=[R2m tMUA_max.R2(s)];
SZ=[SZ tMUA_max.R2(s) tMUA_max.rfs(s)];
    
model='classicRF';
s = strcmp(tMUA_max.Model,model);
C=[C tMUA_max.X(s)./668.745 tMUA_max.Y(s)./668.745 tMUA_max.rfs(s)./2]; 

model='css_ephys_cv1';
s = strcmp(tLFP_max.Model,model);
sig=unique(tLFP_max.SigType);
lfp_order = [3 1 2 5 4];
for i=lfp_order
    b=strcmp(tLFP_max.SigType,sig{i});
    C=[C tLFP_max.R2(s & b) tLFP_max.X(s & b) tLFP_max.Y(s & b) tLFP_max.rfs(s & b)];
    R2m=[R2m tLFP_max.R2(s & b)];
    SZ=[SZ tLFP_max.R2(s & b) tLFP_max.rfs(s & b)];
end

s= sum(R2m>rth,2)==size(R2m,2);
distRF = [...
    sqrt(((C(s,2)-C(s,5)).^2) + ((C(s,3)-C(s,6)).^2)) ...
    sqrt(((C(s,2)-C(s,9)).^2) + ((C(s,3)-C(s,10)).^2)) ...
    sqrt(((C(s,2)-C(s,13)).^2) + ((C(s,3)-C(s,14)).^2)) ...
    sqrt(((C(s,2)-C(s,17)).^2) + ((C(s,3)-C(s,18)).^2)) ...
    sqrt(((C(s,2)-C(s,21)).^2) + ((C(s,3)-C(s,22)).^2)) ];

distSZ = [...
    C(s,4)-C(s,7) ...
    C(s,4)-C(s,11) ...
    C(s,4)-C(s,15) ...
    C(s,4)-C(s,19) ...
    C(s,4)-C(s,23) ];

distSZ = [...
    C(s,4)-C(s,7) ...
    C(s,4)-C(s,11) ...
    C(s,4)-C(s,15) ...
    C(s,4)-C(s,19) ...
    C(s,4)-C(s,23) ];

normSz = [...
    C(s,7)./C(s,4) ...
    C(s,11)./C(s,4) ...
    C(s,15)./C(s,4) ...
    C(s,19)./C(s,4) ...
    C(s,23)./C(s,4) ];

%% MUA model comparison ===================================================
m=unique(tMUA_max.Model);
R2=[];
for i=1:length(m)
    R2 = [R2 tMUA_max.R2(strcmp(tMUA_max.Model,m{i}))];
end

f=figure;
set(f,'Position',[100 100 1200 1200]);
for row=1:4
    for column=1:4
        subplot(4,4,((row-1)*4)+column); hold on;
        plot([0 100],[0 100],'k');
        scatter(R2(:,row+1), R2(:,column+1),'Marker','.',...
            'MarkerEdgeColor',[.3 .3 .3]);
        set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
        xlabel(m{row+1},'interpreter','none'); 
        ylabel(m{column+1},'interpreter','none');
        title('MUA');
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
    for i=1:length(m)
        R2 = [R2 tLFP_max.R2(...
            strcmp(tLFP_max.Model,m{i}) & ...
            strcmp(tLFP_max.SigType,sig{fb}))];
    end
    
    subplot(length(sig),4,spn); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(:,3), R2(:,1),'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    xlabel(m{3},'interpreter','none');ylabel(m{1},'interpreter','none')
    title(sig{fb})
    spn=spn+1;
    
    subplot(length(sig),4,spn); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(:,3), R2(:,2),'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    xlabel(m{3},'interpreter','none');ylabel(m{2},'interpreter','none')
    title(sig{fb})
    spn=spn+1;
    
    subplot(length(sig),4,spn); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(:,1), R2(:,2),'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    xlabel(m{1},'interpreter','none');ylabel(m{2},'interpreter','none')
    title(sig{fb})
    spn=spn+1;
    
    subplot(length(sig),4,spn); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(:,3), R2(:,4),'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    xlabel(m{3},'interpreter','none');ylabel(m{4},'interpreter','none')
    title(sig{fb})
    spn=spn+1;
end

%% pRF size for different ephys signals ===================================
% SZ is [ MUA_R2(1) MUA_RFS(2) 
%         THETA_R2(3) THETA_RFS(4) 
%         ALPHA_R2(5) ALPHA_RFS(6) 
%         BETA_R2(7) BETA_RFS8) 
%         LGAM_R2(9) LG_RFS(10) 
%         HGAM_R2(11) HGAM_RFS(12)]

LAB=['MUA';sig(lfp_order)];

f=figure; set(f,'Position',[100 100 1300 1200]);
r2th=25;

c=0;d=0;
for ref=1:2:12
    c=c+1;
    SS(c).nSZ =[];
    for fb=1:2:12
        d=d+1;
        s=(SZ(:,ref)>r2th & SZ(:,fb)>r2th);
        
        subplot(6,6,d); hold on; plot([0 100],[0 100],'k');
        
        scatter(SZ(s,ref+1),SZ(s,fb+1),'MarkerEdgeColor','k')
        xlabel(LAB{(ref+1)/2});ylabel(LAB{(fb+1)/2});
        title('pRF size')
        
        set(gca,'xlim',[0 10],'ylim',[0 10]);
        
        SS(c).nSZ = [SS(c).nSZ ; ...
            median(diff(SZ(s,[ref+1 fb+1]),1,2)) ...
            std(diff(SZ(s,[ref+1 fb+1]),1,2))./sqrt(sum(s)) ...
            mean(  SZ(s,fb+1)./SZ(s,ref+1) ) ...
            std(  SZ(s,fb+1)./SZ(s,ref+1) )./sqrt(sum(s))];
    end
end

%% What's specific about the good DoG fits=================================



%% Pattern space Size


