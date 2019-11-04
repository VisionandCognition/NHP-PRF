% ck_MakePlots
% Takes the pre-processed fitting result tables and creates comparison plots

%% Paths ==================================================================
BaseFld = pwd;
ResFld = ...
    '/Users/chris/Documents/MRI_ANALYSIS/NHP-PRF/FitResults/MultiModal';
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
end


%% scatter plots & differences -----
figure;

s_R2 = T(7).mod.R2>0;
roi={'V1','V2_merged','V3_merged','V4_merged','MT','MST','TEO','LIP_merged'};


subplot(1,3,1); hold on;
plot([0 100],[0 100],'k');
scatter(T(7).mod.R2(s_R2),T(2).mod.R2(s_R2),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(7).mod.R2(s_R2 & T(7).mod.(roi{r})),...
        T(2).mod.R2(s_R2 & T(2).mod.(roi{r})),'Marker','.');
end
title('linear vs css'); xlabel('linear');ylabel('css');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

subplot(1,3,2); hold on;
plot([0 100],[0 100],'k'); 
scatter(T(7).mod.R2(s_R2),T(4).mod.R2(s_R2),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(7).mod.R2(s_R2 & T(7).mod.(roi{r})),...
        T(4).mod.R2(s_R2 & T(4).mod.(roi{r})),'Marker','.');
end
title('linear vs dog'); xlabel('linear');ylabel('dog');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

% subplot(2,2,3); hold on;
% plot([0 100],[0 100],'k'); 
% scatter(T(7).mod.R2(s_R2),T(8).mod.R2(s_R2),'Marker','.',...
%     'MarkerEdgeColor',[.3 .3 .3]);
% for r=1:length(roi)
%     scatter(T(7).mod.R2(s_R2 & T(7).mod.(roi{r})),...
%         T(8).mod.R2(s_R2 & T(8).mod.(roi{r})),'Marker','.');
% end
% title('linear vs linear_neg'); xlabel('linear');ylabel('linear_neg');
% set(gca, 'Box','on', 'xlim', [0 100], 'ylim',[0 100]);

subplot(1,3,3); hold on;
plot([0 100],[0 100],'k'); 
scatter(T(2).mod.R2(s_R2),T(4).mod.R2(s_R2),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(2).mod.R2(s_R2 & T(2).mod.(roi{r})),...
        T(4).mod.R2(s_R2 & T(4).mod.(roi{r})),'Marker','.');
end
title('css vs dog'); xlabel('css');ylabel('dog');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

%% diff distrutions plots -----
%s_R2 = T(7).mod.R2>0;
%roi={'V1','V2_merged','V3_merged','V4_merged','MT','TEO'};

figure;
%subplot(2,3,4); hold on;

diffmat{1}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(2).mod.R2(s_R2 & T(2).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(2).mod.(roi{r}));
    
    m = mean(T(2).mod.R2(s_R2 & T(2).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})));
    sd = std(T(2).mod.R2(s_R2 & T(2).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(2).mod.(roi{r})));
    diffmat{1} = [diffmat{1}; m sd se];
    
    %bar(x,f);
end

%subplot(2,2,2); hold on;

diffmat{2}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(4).mod.R2(s_R2 & T(4).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(4).mod.(roi{r}));
    
    m = mean(T(4).mod.R2(s_R2 & T(4).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})));
    sd = std(T(4).mod.R2(s_R2 & T(4).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(4).mod.(roi{r})));
    diffmat{2} = [diffmat{2}; m sd se];
    
    %bar(x,f);
end

%subplot(2,2,3); hold on;

diffmat{3}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(8).mod.R2(s_R2 & T(8).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(8).mod.(roi{r}));
    
    m = mean(T(8).mod.R2(s_R2 & T(8).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})));
    sd = std(T(8).mod.R2(s_R2 & T(8).mod.(roi{r}))-...
        T(7).mod.R2(s_R2 & T(7).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(8).mod.(roi{r})));
    diffmat{3} = [diffmat{3}; m sd se];
    
    %bar(x,f);
end

%subplot(2,2,4); hold on;

diffmat{3}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(4).mod.R2(s_R2 & T(4).mod.(roi{r}))-...
        T(2).mod.R2(s_R2 & T(2).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(2).mod.(roi{r}));
    
    m = mean(T(4).mod.R2(s_R2 & T(4).mod.(roi{r}))-...
        T(2).mod.R2(s_R2 & T(2).mod.(roi{r})));
    sd = std(T(4).mod.R2(s_R2 & T(4).mod.(roi{r}))-...
        T(2).mod.R2(s_R2 & T(2).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(2).mod.(roi{r})));
    diffmat{3} = [diffmat{3}; m sd se];
    
    %bar(x,f);
end

figure;
for cc = 1:length(diffmat)
    subplot(1,3,cc); hold on;
    bar(1:length(diffmat{cc}),diffmat{cc}(:,1));
    errorbar(1:length(diffmat{cc}),diffmat{cc}(:,1),diffmat{cc}(:,3),'Linestyle','none')
end



%% scatter plot HRF & differences -----
figure;

s_R2 = T(7).mod.R2>0;
roi={'V1','V2_merged','V3_merged','V4_merged','MT','MST','TEO','LIP_merged'};

hold on;
plot([0 100],[0 100],'k');
scatter(T(1).mod.R2(s_R2),T(2).mod.R2(s_R2),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(1).mod.R2(s_R2 & T(1).mod.(roi{r})),...
        T(2).mod.R2(s_R2 & T(2).mod.(roi{r})),'Marker','.');
end
title('Canonical vs monkey HRF'); xlabel('Monkey HRF R2');ylabel('Canonical HRF R2');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

diffmat2{1}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(2).mod.R2(s_R2 & T(2).mod.(roi{r}))-...
        T(1).mod.R2(s_R2 & T(1).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(2).mod.(roi{r}));
    
    m = mean(T(2).mod.R2(s_R2 & T(2).mod.(roi{r}))-...
        T(1).mod.R2(s_R2 & T(1).mod.(roi{r})));
    sd = std(T(2).mod.R2(s_R2 & T(2).mod.(roi{r}))-...
        T(1).mod.R2(s_R2 & T(1).mod.(roi{r})));
    se = sd ./ sqrt(sum(s_R2 & T(2).mod.(roi{r})));
    diffmat2{1} = [diffmat2{1}; m sd se];
    
    %bar(x,f);
end

figure; hold on
bar(1:length(diffmat2{1}),diffmat2{1}(:,1));
errorbar(1:length(diffmat2{1}),diffmat2{1}(:,1),diffmat2{1}(:,3),'Linestyle','none')

%% rf size =====
figure;
s_R2 = T(2).mod.R2>5;
roi={'V1','V2_merged','V3_merged','V4_merged','MT','MST','TEO','LIP_merged'};

hold on;
plot([0 100],[0 100],'k');
scatter(T(1).mod.rfs(s_R2),T(2).mod.rfs(s_R2),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(1).mod.rfs(s_R2 & T(1).mod.(roi{r})),...
        T(2).mod.rfs(s_R2 & T(2).mod.(roi{r})),'Marker','.');
end
title('Canonical vs monkey HRF'); xlabel('Canonical HRF sz');ylabel('Monkey HRF sz');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

diffmat2{1}=[];
for r=1:length(roi)
    [n,x] = hist(...
        T(2).mod.rfs(s_R2 & T(2).mod.(roi{r}))-...
        T(1).mod.rfs(s_R2 & T(1).mod.(roi{r})),...
        100);
    f = n./sum(s_R2 & T(2).mod.(roi{r}));
    
    dsz = T(2).mod.rfs-T(1).mod.rfs;
    dsz = dsz(s_R2 & T(2).mod.(roi{r}));
    dsz = dsz(isfinite(dsz));
    
    m = mean(dsz);
    sd = std(dsz);
    se = sd ./ sqrt(length(dsz));
    diffmat2{1} = [diffmat2{1}; m sd se];
    
    %bar(x,f);
end

figure; hold on
bar(1:length(diffmat2{1}),diffmat2{1}(:,1));
errorbar(1:length(diffmat2{1}),diffmat2{1}(:,1),diffmat2{1}(:,3),'Linestyle','none')











%% ECC vs Size plot s -----
figure;

s_R2 = T(2).mod.R2>4;
roi={'V1','V2_merged','V3_merged','V4_merged','MT','MST','TEO','LIP_merged'};

hold on;
scatter(T(2).mod.ecc(s_R2),T(2).mod.rfs(s_R2),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
for r=1:length(roi)
    scatter(T(2).mod.ecc(s_R2 & T(2).mod.(roi{r})),...
        T(2).mod.rfs(s_R2 & T(2).mod.(roi{r})),'Marker','.');
end
title('ecc vs rfs'); xlabel('ecc');ylabel('rfs');
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

%% Ephys VFC ----
% MUA
s=tMUA_max.R2>50;
szm=tMUA_max.rfs<5;

model='css_ephys_cv1';

figure;

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
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

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
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

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
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

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
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

%% Ephys VFC ----
% MUA
s=tLFP_max.R2>10;

model='css_ephys_cv1';
b = 'lGamma';
figure;

subplot(2,2,1);hold on;
szm=tLFP_max.rfs<5;
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
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

subplot(2,2,3);hold on;
szm=tLFP_max.rfs<5;
m=strcmp(tLFP_max.Monkey,'lick') & strcmp(tLFP_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tLFP_max.Area==4;
for r=unique(tLFP_max.Array)'
    a=tLFP_max.Array==r;
    scatter(tLFP_max.X(s & m & a & v),...
        tLFP_max.Y(s & m & a & v),'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

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
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

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
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);



%% Ephys VFC ----
% LFP gamma
s=tMUA_max.SNR>2;
model='classicRF';

figure;

subplot(2,2,1);hold on;
m=strcmp(tMUA_max.Monkey,'lick') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==1;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v)./668.745,...
        tMUA_max.Y(s & m & a & v)./668.745,'Marker','*' )
    sum(s & m & a & v)
end
%set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

subplot(2,2,3);hold on;
m=strcmp(tMUA_max.Monkey,'lick') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==4;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v)./668.745,...
        tMUA_max.Y(s & m & a & v)./668.745,'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

subplot(2,2,2);hold on;
m=strcmp(tMUA_max.Monkey,'aston') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==1;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v)./668.745,...
        tMUA_max.Y(s & m & a & v)./668.745,'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

subplot(2,2,4);hold on;
m=strcmp(tMUA_max.Monkey,'aston') & strcmp(tMUA_max.Model,model);
plot([-10 20],[0 0],'k');
plot([0 0],[-30 10],'k');
v=tMUA_max.Area==4;
for r=unique(tMUA_max.Array)'
    a=tMUA_max.Array==r;
    scatter(tMUA_max.X(s & m & a & v)./668.745,...
        tMUA_max.Y(s & m & a & v)./668.745,'Marker','*' )

end
set(gca, 'Box','off', 'xlim', [-10 20], 'ylim',[-30 10]);

%% Ephys location difference & size difference  ---
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
sig=unique(tLFP_max.SigType)
for i=[2 5 4]
    b=strcmp(tLFP_max.SigType,sig{i});
    C=[C tLFP_max.R2(s & b) tLFP_max.X(s & b) tLFP_max.Y(s & b) tLFP_max.rfs(s & b)];
    R2m=[R2m tLFP_max.R2(s & b)];
    SZ=[SZ tLFP_max.R2(s & b) tLFP_max.rfs(s & b)];
end


s= sum(R2m>50,2)==size(R2m,2);
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


%% MUA model comparison ---
m=unique(tMUA_max.Model);
R2=[];
for i=1:length(m)
    R2 = [R2 tMUA_max.R2(strcmp(tMUA_max.Model,m{i}))];
end

figure;
subplot(1,3,1); hold on;
plot([0 100],[0 100],'k');
scatter(R2(:,4), R2(:,2),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

subplot(1,3,2); hold on;
plot([0 100],[0 100],'k');
scatter(R2(:,4), R2(:,3),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

subplot(1,3,3); hold on;
plot([0 100],[0 100],'k');
scatter(R2(:,2), R2(:,3),'Marker','.',...
    'MarkerEdgeColor',[.3 .3 .3]);
set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);

%% LFP model comparison ---
l=unique(tLFP_max.SigType);
for fb=1:length(l)
    
    m=unique(tLFP_max.Model);
    R2=[];
    for i=1:length(m)
        R2 = [R2 tLFP_max.R2(...
            strcmp(tLFP_max.Model,m{i}) & ...
            strcmp(tLFP_max.SigType,l{fb}))];
    end
    
    figure;
    subplot(1,3,1); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(:,4), R2(:,2),'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    
    subplot(1,3,2); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(:,4), R2(:,3),'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
    
    subplot(1,3,3); hold on;
    plot([0 100],[0 100],'k');
    scatter(R2(:,2), R2(:,3),'Marker','.',...
        'MarkerEdgeColor',[.3 .3 .3]);
    set(gca, 'Box','off', 'xlim', [0 100], 'ylim',[0 100]);
end

%%
r2th=50;
c=0;
for ref=1:2:7
    c=c+1;
    figure;
    d=0;
    SS(c).nSZ =[];
    for fb=1:2:9
        d=d+1;
        s=(SZ(:,ref)>r2th & SZ(:,fb)>r2th & SZ(:,ref+1)<10 & SZ(:,fb+1)<10);
        subplot(1,5,d); hold on; plot([0 10],[0 10],'k');
        scatter(SZ(s,ref+1),SZ(s,fb+1))
        SS(c).nSZ = [SS(c).nSZ ; ...
            median(diff(SZ(s,[ref+1 fb+1]),1,2)) std(diff(SZ(s,[ref+1 fb+1]),1,2))./sqrt(sum(s)) ...
            mean(  SZ(s,fb+1)./SZ(s,ref+1) ) std(  SZ(s,fb+1)./SZ(s,ref+1) )./sqrt(sum(s))];
    end
end



%%
s=(SZ(:,1)>20 & SZ(:,3)>20 & SZ(:,2)<10 & SZ(:,6)<10);
nSZ = [nSZ ; ...
    mean(diff(SZ(s,[2 6]),1,2)) std(diff(SZ(s,[2 6]),1,2))./sqrt(sum(s)) ...
    mean(  SZ(s,6)./SZ(s,2) ) std(  SZ(s,6)./SZ(s,2) )./sqrt(sum(s))];

s=(SZ(:,1)>20 & SZ(:,5)>20 & SZ(:,2)<10 & SZ(:,6)<10);
nSZ = [nSZ ; ...
    mean(diff(SZ(s,[2 6]),1,2)) std(diff(SZ(s,[2 6]),1,2))./sqrt(sum(s)) ...
    mean(  SZ(s,6)./SZ(s,2) ) std(  SZ(s,6)./SZ(s,2) )./sqrt(sum(s))];

s=(SZ(:,1)>20 & SZ(:,7)>20 & SZ(:,2)<10 & SZ(:,6)<10);
nSZ = [nSZ ; ...
    mean(diff(SZ(s,[2 6]),1,2)) std(diff(SZ(s,[2 6]),1,2))./sqrt(sum(s)) ...
    mean(  SZ(s,6)./SZ(s,2) ) std(  SZ(s,6)./SZ(s,2) )./sqrt(sum(s))];

s=(SZ(:,1)>20 & SZ(:,9)>20 & SZ(:,2)<10 & SZ(:,6)<10);
nSZ = [nSZ ; ...
    mean(diff(SZ(s,[2 6]),1,2)) std(diff(SZ(s,[2 6]),1,2))./sqrt(sum(s)) ...
    mean(  SZ(s,6)./SZ(s,2) ) std(  SZ(s,6)./SZ(s,2) )./sqrt(sum(s))];