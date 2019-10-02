%%
load pRF_estimates_MRI_inFunc

%%
VOX = []; 
ROImap={...
    'V1','V2','V3','V3A','V4','MT','MST',...
    'LIP','VIP','TEO','TPO','Area5','Area7',...
    'V1_elec','V4_elec'};

for m=1:2
    nv=numel(MRI(m).V1.img);
    VOX(m).POL = reshape(MRI(m).prf_pol.img,nv,1);
    VOX(m).ECC = reshape(MRI(m).prf_ecc.img,nv,1);
    VOX(m).R2 = reshape(MRI(m).prf_R2.img,nv,1);
    VOX(m).SZ = reshape(MRI(m).prf_sz.img,nv,1);
    VOX(m).SD = reshape(MRI(m).prf_sd.img,nv,1);
    VOX(m).ROI = logical( [...
        reshape(MRI(m).V1.img,nv,1) ...
        reshape(MRI(m).V2.img,nv,1) ...
        reshape(MRI(m).V3.img,nv,1) ...
        reshape(MRI(m).V3A.img,nv,1) ...
        reshape(MRI(m).V4.img,nv,1) ...
        reshape(MRI(m).MT.img,nv,1) ...
        reshape(MRI(m).MST.img,nv,1) ...
        reshape(MRI(m).LIP.img,nv,1) ...
        reshape(MRI(m).VIP.img,nv,1) ...
        reshape(MRI(m).TEO.img,nv,1) ...
        reshape(MRI(m).TPO.img,nv,1) ...
        reshape(MRI(m).Area5.img,nv,1) ...
        reshape(MRI(m).Area7.img,nv,1) ...
        reshape(MRI(m).V1Electr.img,nv,1) ...
        reshape(MRI(m).V4Electr.img,nv,1) ] );
    VOX(m).brainmask=reshape(MRI(m).brainmask.img,nv,1);
end    


% %%  
% %R2 distributions
% for a=1:15
%     r2 = [VOX(1).R2(VOX(1).ROI(:,a));VOX(2).R2(VOX(2).ROI(:,a))];
%     [N,X]=hist(r2,100);    
%     
%     figure;
%     bar(X,N)
%     title(ROImap{a})
%     xlabel('R2'); ylabel('Proportion of voxels');
% end

%% V1 distributions of SZ, ECC and POL
m=1; a=1; th=5;
max_sz = 30;
max_ecc = 10;

ecc = [VOX(m).ECC(VOX(m).ROI(:,a) & VOX(m).R2>th & ...
    VOX(m).ECC<max_ecc & VOX(m).SZ<max_sz)];
sz = [VOX(m).SZ(VOX(m).ROI(:,a) & VOX(m).R2>th & ...
    VOX(m).ECC<max_ecc & VOX(m).SZ<max_sz)];
sd = [VOX(m).SD(VOX(m).ROI(:,a) & VOX(m).R2>th & ...
    VOX(m).ECC<max_ecc & VOX(m).SD<max_sz)];
pol = [VOX(m).POL(VOX(m).ROI(:,a) & VOX(m).R2>th & ...
    VOX(m).ECC<max_ecc & VOX(m).SZ<max_sz)];

figure;
subplot(1,3,1)
hist(ecc,100)
subplot(1,3,2)
hist(sz,40)
subplot(1,3,3)
hist(pol,100)

%% ECC distributions only
figure;
m=2;
th=8;
[Y,X]=hist(VOX(m).ECC(VOX(m).R2>th),0:0.5:100);
bar(X,Y./sum(VOX(m).R2>th))

%% VFC all voxels
th=5;
for m=1:2
    [X,Y]=pol2cart(VOX(m).POL.*(pi/180), VOX(m).ECC);
    R2=VOX(m).R2;
    SD=VOX(m).SD;
    
    figure(m);
    
    inc=R2>th;
    
    subplot(1,2,1);hold on;
    plot(X(inc),Y(inc),'*','MarkerSize',8);
    set(gca,'xlim',[-15 15],'ylim',[-15 15],'FontSize',14);
    title(['M' num2str(m) ' ALL'],'FontSize', 18);
    viscircles([0 0], 8,'LineWidth',2,'Color','k');
    
    subplot(1,2,2);hold on;
    viscircles([X(inc),Y(inc)],SD(inc)./2,'LineWidth',1);
    set(gca,'xlim',[-15 15],'ylim',[-15 15],'FontSize',14);
    title(['M' num2str(m) ' ALL'],'FontSize', 18);
    viscircles([0 0], 8,'LineWidth',2,'Color','k');
      
    for j=1
        subplot(1,2,j);
        text(13.5,-13.5,['R2 > ' num2str(th)],...
            'HorizontalALignment','right','FontSize',14)
        text(13.5,13.5,['n = ' num2str(sum(inc))],...
            'HorizontalALignment','right','FontSize',14)
    end
end

%% R2 distribution
figure
for m=1:2
    subplot(1,2,m)
    R2=VOX(m).R2(logical(VOX(m).brainmask));
    R2=R2(R2>0);
    [Y,X]=hist(R2,1:1:100);
    bar(X,Y/length(R2));
end

%%
figure;
q=1;
for a=1:2
    for m=1:2
        subplot(2,2,q)
        R2=VOX(m).R2(logical(VOX(m).ROI(:,a)));
        R2=R2(R2>0);
        [Y,X]=hist(R2,1:1:100);
        bar(X,Y/length(R2));
        q=q+1;
        
        mean(R2)
        max(R2)
    end
end




%% VFC all V1 & V4 voxels
th=5;
figure;
for m=1:2
    inc=VOX(m).ROI(:,1);
    [X,Y]=pol2cart(VOX(m).POL(inc).*(pi/180), VOX(m).ECC(inc));
    R2=VOX(m).R2(inc);
    SD=VOX(m).SD(inc);
    inc=R2>th;
    
    subplot(2,2,m);hold on;
    plot(X(inc),Y(inc),'*','MarkerSize',8);
    set(gca,'xlim',[-15 15],'ylim',[-15 15],'FontSize',14);
    title(['M' num2str(m) ' V1'],'FontSize', 18);
    viscircles([0 0], 8,'LineWidth',2,'Color','k');
    
    inc=VOX(m).ROI(:,5);
    [X,Y]=pol2cart(VOX(m).POL(inc).*(pi/180), VOX(m).ECC(inc));
    R2=VOX(m).R2(inc);
    SD=VOX(m).SD(inc);
    inc=R2>th;
    
    subplot(2,2,m+2);hold on;
    plot(X(inc),Y(inc),'*','MarkerSize',8);
    %viscircles([X(inc),Y(inc)],SD(inc)./2,'LineWidth',1);
    set(gca,'xlim',[-15 15],'ylim',[-15 15],'FontSize',14);
    title(['M' num2str(m) ' V4'],'FontSize', 18);
    viscircles([0 0], 8,'LineWidth',2,'Color','k');
      
%     for j=1
%         subplot(2,2,j);
%         text(13.5,-13.5,['R2 > ' num2str(th)],...
%             'HorizontalALignment','right','FontSize',14)
%         text(13.5,13.5,['n = ' num2str(sum(inc))],...
%             'HorizontalALignment','right','FontSize',14)
%     end
end

%% ECC vs Size
th=5;
ecc_range=[-2 8];
f1=figure;
f2=figure;hold on
for a=1:13
    figure(f1);
    subplot(4,4,a);hold on;
    m=1;
    ecc1 = [VOX(m).ECC(VOX(m).ROI(:,a) & VOX(m).R2>th)];...
    sz1 = [VOX(m).SZ(VOX(m).ROI(:,a) & VOX(m).R2>th)];...
    sd1 = [VOX(m).SD(VOX(m).ROI(:,a) & VOX(m).R2>th)];...
    m=2;
    ecc2 = [VOX(m).ECC(VOX(m).ROI(:,a) & VOX(m).R2>th)];...
    sz2 = [VOX(m).SZ(VOX(m).ROI(:,a) & VOX(m).R2>th)];...
    sd2 = [VOX(m).SD(VOX(m).ROI(:,a) & VOX(m).R2>th)];...

    ecc = [ecc1;ecc2];
    sd = [sd1;sd2];
    sz = [sz1; sz2]; 
%     % bin it
%     ec_sz=[];
%     for i=0.5:0.5:ecc_max
%         idx = (ecc > i-0.5 & ecc < i);
%         ec_sz=[ec_sz; ...
%             i-0.25 nanmean(sd(idx)) median(sd(idx)) nanstd(sd(idx)) sum(idx)];
%     end
    
    %figure;
    scatter(ecc,sd)
    title(ROImap{a});
    ylabel('pRF Size (sigma)'); xlabel('Eccentricity (deg)');
    if a==5
        ecc_range=ecc_range+4;
    end
    p=polyfit(ecc(ecc>ecc_range(1) & ecc<ecc_range(2)),sd(ecc>ecc_range(1) & ecc<ecc_range(2)),1);
    x=[0 24];y=polyval(p,x);
    plot(x,y, 'LineWidth',3)
    
    set(gca,'Ylim',[-.5 5])
    
%     hold on;
%     [p1,p2]=polyfit(ec_sz(1:10,1),ec_sz(1:10,2),1);
%     xf=[0 ecc_max];
%     f=polyval(p1,xf);
%     plot(xf,f);
%     
%     legend(ROImap)

    figure(f2);
    plot(x,y, 'LineWidth',3)
    
    if a==1
        poly_MRI_V1 = p;
    elseif a==5
        poly_MRI_V4 = p;
    end
    
end
figure(f2);
legend(ROImap)
