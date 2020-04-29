addpath(genpath('../../FitResults/MultiModal'));

load pRF_estimates_MRI_inFunc
load pRF_estimates_ephys

%%
MUA=[RetMap(1).table_mua;RetMap(2).table_mua];
LFP=[RetMap(1).table_lfp;RetMap(2).table_lfp];

VOX = []; 
ROImap={...
    'V1','V4',...
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
        reshape(MRI(m).V4.img,nv,1) ...
        reshape(MRI(m).V1Electr.img,nv,1) ...
        reshape(MRI(m).V4Electr.img,nv,1) ] );
    VOX(m).brainmask=reshape(MRI(m).brainmask.img,nv,1);
    
    [VOX(m).X, VOX(m).Y] = pol2cart(...
        deg2rad(VOX(m).POL), VOX(m).ECC );
end    

%%
th_mri=5;
th_ephys=.6;

extend=5;
x=1:extend;
y=-1:-1:-extend;
mS_mri_V1=zeros(extend);
mS_mri_V4=zeros(extend);

mS_mua_V1=zeros(extend);
mS_mua_V4=zeros(extend);

for xi=x
    for yi=y
        
        X=[VOX(1).X;VOX(2).X];
        Y=[VOX(1).Y;VOX(2).Y];
        
        E=[VOX(1).ECC;VOX(2).ECC];
        S=[VOX(1).SD;VOX(2).SD];
        
        ROI=[VOX(1).ROI;VOX(2).ROI];

        R2=[VOX(1).R2;VOX(2).R2];
        
        idx = (X>xi-1 & X<xi & Y>yi & Y<yi+1);
        idx2 = R2>th_mri;
        
        idx = logical(idx.*idx2);
        
        coll(abs(yi),xi).E = E(idx);
        coll(abs(yi),xi).S = S(idx);
        coll(abs(yi),xi).X = X(idx);
        coll(abs(yi),xi).Y = Y(idx);
        coll(abs(yi),xi).ROI = ROI(idx,:);
        
        mS_mri_V1(abs(yi),xi) = nanmean(coll(abs(yi),xi).S(coll(abs(yi),xi).ROI(:,1),:));
        mS_mri_V4(abs(yi),xi) = nanmean(coll(abs(yi),xi).S(coll(abs(yi),xi).ROI(:,2),:));
        
        mS_mri_V1e(abs(yi),xi) = nanmean(coll(abs(yi),xi).S(coll(abs(yi),xi).ROI(:,3),:));
        mS_mri_V4e(abs(yi),xi) = nanmean(coll(abs(yi),xi).S(coll(abs(yi),xi).ROI(:,4),:));
        
        % =======
        
        X=MUA.prf_x;
        Y=MUA.prf_y;
        
        E=MUA.prf_ecc;
        S=MUA.prf_sd;
        R2=MUA.prf_r2;
        
        idx = (X>xi-1 & X<xi & Y>yi & Y<yi+1);
        idx2 = R2>th_ephys;
        
        idx = logical(idx.*idx2);
        
        coll_mua(abs(yi),xi).E = E(idx);
        coll_mua(abs(yi),xi).S = S(idx);
        coll_mua(abs(yi),xi).X = X(idx);
        coll_mua(abs(yi),xi).Y = Y(idx);
        
        mS_mua_V1(abs(yi),xi) = nanmean(coll_mua(abs(yi),xi).S(MUA.Area(idx)==1));
        mS_mua_V4(abs(yi),xi) = nanmean(coll_mua(abs(yi),xi).S(MUA.Area(idx)==4));
        
        % =======
        
        X=LFP.alpha_prf_x;
        Y=LFP.alpha_prf_y;
        
        E=LFP.alpha_prf_ecc;
        S=LFP.alpha_prf_sd;
        R2=LFP.alpha_prf_r2;
        
        idx = (X>xi-1 & X<xi & Y>yi & Y<yi+1);
        idx2 = R2>th_ephys;
        
        idx = logical(idx.*idx2);
        
        coll_lfpa(abs(yi),xi).E = E(idx);
        coll_lfpa(abs(yi),xi).S = S(idx);
        coll_lfpa(abs(yi),xi).X = X(idx);
        coll_lfpa(abs(yi),xi).Y = Y(idx);
        
        mS_lfpa_V1(abs(yi),xi) = nanmean(coll_lfpa(abs(yi),xi).S(LFP.Area(idx)==1));
        mS_lfpa_V4(abs(yi),xi) = nanmean(coll_lfpa(abs(yi),xi).S(LFP.Area(idx)==4));

        % =======
        
        X=LFP.beta_prf_x;
        Y=LFP.beta_prf_y;
        
        E=LFP.beta_prf_ecc;
        S=LFP.beta_prf_sd;
        R2=LFP.beta_prf_r2;

        idx = (X>xi-1 & X<xi & Y>yi & Y<yi+1);
        idx2 = R2>th_ephys;
        
        idx = logical(idx.*idx2);
        
        coll_lfpb(abs(yi),xi).E = E(idx);
        coll_lfpb(abs(yi),xi).S = S(idx);
        coll_lfpb(abs(yi),xi).X = X(idx);
        coll_lfpb(abs(yi),xi).Y = Y(idx);
        
        mS_lfpb_V1(abs(yi),xi) = nanmean(coll_lfpb(abs(yi),xi).S(LFP.Area(idx)==1));
        mS_lfpb_V4(abs(yi),xi) = nanmean(coll_lfpb(abs(yi),xi).S(LFP.Area(idx)==4));
        
        % =======
        
        X=LFP.lowgamma_prf_x;
        Y=LFP.lowgamma_prf_y;
        
        E=LFP.lowgamma_prf_ecc;
        S=LFP.lowgamma_prf_sd;
        R2=LFP.lowgamma_prf_r2;
        
        idx = (X>xi-1 & X<xi & Y>yi & Y<yi+1);
        idx2 = R2>th_ephys;
        
        idx = logical(idx.*idx2);
        
        coll_lfpg1(abs(yi),xi).E = E(idx);
        coll_lfpg1(abs(yi),xi).S = S(idx);
        coll_lfpg1(abs(yi),xi).X = X(idx);
        coll_lfpg1(abs(yi),xi).Y = Y(idx);
        
        mS_lfpg1_V1(abs(yi),xi) = nanmean(coll_lfpg1(abs(yi),xi).S(LFP.Area(idx)==1));
        mS_lfpg1_V4(abs(yi),xi) = nanmean(coll_lfpg1(abs(yi),xi).S(LFP.Area(idx)==4));
        
        % =======
        
        X=LFP.highgamma_prf_x;
        Y=LFP.highgamma_prf_y;
        
        E=LFP.highgamma_prf_ecc;
        S=LFP.highgamma_prf_sd;
        R2=LFP.highgamma_prf_r2;
        
        idx = (X>xi-1 & X<xi & Y>yi & Y<yi+1);
        idx2 = R2>th_ephys;
        
        idx = logical(idx.*idx2);
        
        coll_lfpg2(abs(yi),xi).E = E(idx);
        coll_lfpg2(abs(yi),xi).S = S(idx);
        coll_lfpg2(abs(yi),xi).X = X(idx);
        coll_lfpg2(abs(yi),xi).Y = Y(idx);
        
        mS_lfpg2_V1(abs(yi),xi) = nanmean(coll_lfpg2(abs(yi),xi).S(LFP.Area(idx)==1));
        mS_lfpg2_V4(abs(yi),xi) = nanmean(coll_lfpg2(abs(yi),xi).S(LFP.Area(idx)==4));
    end
end

%%
figure;
subplot(2,2,1)
imagesc(mS_mri_V1,[0.5 2]);
subplot(2,2,2)
imagesc(mS_mri_V4,[0.5 2]);
subplot(2,2,3)
imagesc(mS_mri_V1e,[0.5 2]);
subplot(2,2,4)
imagesc(mS_mri_V4e,[0.5 2]);

%% 
figure;
subplot(5,2,1)
imagesc(mS_mua_V1,[0.5 3]);
subplot(5,2,2)
imagesc(mS_mua_V4,[0.50 3]);

subplot(5,2,3)
imagesc(mS_lfpa_V1,[0.50 3]);
subplot(5,2,4)
imagesc(mS_lfpa_V4,[0.50 3]);

subplot(5,2,5)
imagesc(mS_lfpb_V1,[0.50 3]);
subplot(5,2,6)
imagesc(mS_lfpb_V4,[0.50 3]);

subplot(5,2,7)
imagesc(mS_lfpg1_V1,[0.50 3]);
subplot(5,2,8)
imagesc(mS_lfpg1_V4,[0.50 3]);

subplot(5,2,9)
imagesc(mS_lfpg2_V1,[0.50 3]);
subplot(5,2,10)
imagesc(mS_lfpg2_V4,[0.50 3]);

%%
ccV1=[];
[R,p] = corrcoef(mS_mri_V1(:) , mS_mua_V1(:), 'rows','complete'); ccV1=[ccV1;R(2)];
[R,p] = corrcoef(mS_mri_V1(:) , mS_lfpa_V1(:), 'rows','complete'); ccV1=[ccV1;R(2)];
[R,p] = corrcoef(mS_mri_V1(:) , mS_lfpb_V1(:), 'rows','complete'); ccV1=[ccV1;R(2)];
[R,p] = corrcoef(mS_mri_V1(:) , mS_lfpg1_V1(:), 'rows','complete'); ccV1=[ccV1;R(2)];
[R,p] = corrcoef(mS_mri_V1(:) , mS_lfpg2_V1(:), 'rows','complete'); ccV1=[ccV1;R(2)];

ccV1e=[];
[R,p] = corrcoef(mS_mri_V1e(:) , mS_mua_V1(:), 'rows','complete'); ccV1e=[ccV1e;R(2)];
[R,p] = corrcoef(mS_mri_V1e(:) , mS_lfpa_V1(:), 'rows','complete'); ccV1e=[ccV1e;R(2)];
[R,p] = corrcoef(mS_mri_V1e(:) , mS_lfpb_V1(:), 'rows','complete'); ccV1e=[ccV1e;R(2)];
[R,p] = corrcoef(mS_mri_V1e(:) , mS_lfpg1_V1(:), 'rows','complete'); ccV1e=[ccV1e;R(2)];
[R,p] = corrcoef(mS_mri_V1e(:) , mS_lfpg2_V1(:), 'rows','complete'); ccV1e=[ccV1e;R(2)];

ccV4=[];
[R,p] = corrcoef(mS_mri_V4(:) , mS_mua_V4(:), 'rows','complete'); ccV4=[ccV4;R(2)];
[R,p] = corrcoef(mS_mri_V4(:) , mS_lfpa_V4(:), 'rows','complete'); ccV4=[ccV4;R(2)];
[R,p] = corrcoef(mS_mri_V4(:) , mS_lfpb_V4(:), 'rows','complete'); ccV4=[ccV4;R(2)];
[R,p] = corrcoef(mS_mri_V4(:) , mS_lfpg1_V4(:), 'rows','complete'); ccV4=[ccV4;R(2)];
[R,p] = corrcoef(mS_mri_V4(:) , mS_lfpg2_V4(:), 'rows','complete'); ccV4=[ccV4;R(2)];

ccV4e=[];
[R,p] = corrcoef(mS_mri_V4e(:) , mS_mua_V4(:), 'rows','complete'); ccV4e=[ccV4e;R(2)];
[R,p] = corrcoef(mS_mri_V4e(:) , mS_lfpa_V4(:), 'rows','complete'); ccV4e=[ccV4e;R(2)];
[R,p] = corrcoef(mS_mri_V4e(:) , mS_lfpb_V4(:), 'rows','complete'); ccV4e=[ccV4e;R(2)];
[R,p] = corrcoef(mS_mri_V4e(:) , mS_lfpg1_V4(:), 'rows','complete'); ccV4e=[ccV4e;R(2)];
[R,p] = corrcoef(mS_mri_V4e(:) , mS_lfpg2_V4(:), 'rows','complete'); ccV4e=[ccV4e;R(2)];

%%
figure;
subplot(1,2,1)
bar(1:5,ccV1)
subplot(1,2,2)
bar(1:5,ccV4)

%%
figure;
x=[0 20];
subplot(1,2,1); hold on;
y=polyval(poly_MRI_V1,x);
plot(x,y)
y=polyval(poly_MUA_V1,x);
plot(x,y)
y=polyval(poly_LFPa_V1,x);
plot(x,y)
y=polyval(poly_LFPb_V1,x);
plot(x,y)
y=polyval(poly_LFPg1_V1,x);
plot(x,y)
y=polyval(poly_LFPg2_V1,x);
plot(x,y)

subplot(1,2,2); hold on;
y=polyval(poly_MRI_V4,x);
plot(x,y)
y=polyval(poly_MUA_V4,x);
plot(x,y)
y=polyval(poly_LFPa_V4,x);
plot(x,y)
y=polyval(poly_LFPb_V4,x);
plot(x,y)
y=polyval(poly_LFPg1_V4,x);
plot(x,y)
y=polyval(poly_LFPg2_V4,x);
plot(x,y)

rmpath(genpath('../../FitResults/MultiModal'));





