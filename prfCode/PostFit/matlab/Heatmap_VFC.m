% Heatmap_VFC.m
% the colorscheme I used depends on:
% https://nl.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps

load('tMUA_css')
tMUA = tMUA_css; 
clear tMUA_css;

%%
R2TH = 50; % R2 inclusion threshold
s=tMUA.R2>R2TH;

%% Heatmap MUA =============================================================
fhm=figure; set(fhm,'Position',[100 100 1800 1000]);
settings.PixPerDeg = 29.5032;
settings.meshsize = 2000;   
% colormap(inferno) % >> see header for dependency

% Lick V1 ----
m=strcmp(tMUA.Monkey,'lick');
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
imagesc(sumimg./max(img(:)));
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
imagesc(sumimg./max(img(:)));
set(gca,'xlim',xrange_idx,'ylim',...
    yrange_idx,'Color','k','xtick',[],'ytick',[])
colorbar;
subplot(2,4,6); hold on;
plot(1.2*res.xr,[0 0],'w'); plot([0 0],1.2*res.yr,'w');
set(gca,'xlim',xrr,'ylim',yrr,'Color','k')
colorbar; clear('res','img','sumimg');
title('Lick V4 MUA')

% Aston V1 ----
m=strcmp(tMUA.Monkey,'aston');
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
imagesc(sumimg./max(img(:)));
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
imagesc(sumimg./max(img(:)));
set(gca,'xlim',xrange_idx,'ylim',...
    yrange_idx,'Color','k','xtick',[],'ytick',[])
colorbar;
subplot(2,4,8); hold on;
plot(1.2*res.xr,[0 0],'w'); plot([0 0],1.2*res.yr,'w');
set(gca,'xlim',xrr,'ylim',yrr,'Color','k')
colorbar; clear('res','img','sumimg');
title('Aston V4 MUA')