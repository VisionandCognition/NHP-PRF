% ck_Figures_Ephys

%% example traces max MUA fit
M=1;
% max r2
idx=find(RetMap(M).table_mua.prf_r2==max(RetMap(M).table_mua.prf_r2))
% mua
muapred = RetMap(M).MUA(RetMap(M).table_mua.Instance(idx)).pred(:,RetMap(M).table_mua.ChanIndex(idx));
mua = RetMap(M).MUA(RetMap(M).table_mua.Instance(idx)).MUAm(:,RetMap(M).table_mua.ChanIndex(idx));
% lfp (same channel)
for f=1:5
    lfppred{f} = RetMap(M).LFP(RetMap(M).table_mua.Instance(idx),f).pred(:,RetMap(M).table_mua.ChanIndex(idx));
    lfp{f} = RetMap(M).LFP(RetMap(M).table_mua.Instance(idx),f).LFPm(:,RetMap(M).table_mua.ChanIndex(idx));
end

figure; 
subplot(6,1,1);hold on;
plot(muapred);
plot(mua, 'o');

for f=1:5
    subplot(6,1,1+f); hold on;
    plot(lfppred{f});
    plot(lfp{f}, 'o');
end

%% Goodness of Fit MUA
figure;
for M=1:2
    subplot(1,2,M);
    r2=RetMap(M).table_mua.prf_r2;
    [Y,X]=hist(r2,0:0.05:1);
    bar(X,Y/length(r2),'BarWidth',1)
end

%% ECC vs Size MUA
figure;th=0.7;
v1_eccsz=[];v4_eccsz=[];

for M=1:2
    subplot(1,2,M);hold on
    r2=RetMap(M).table_mua.prf_r2(RetMap(M).table_mua.prf_r2>th);
    ecc=RetMap(M).table_mua.prf_ecc(RetMap(M).table_mua.prf_r2>th);
    sd=RetMap(M).table_mua.prf_sd(RetMap(M).table_mua.prf_r2>th);
    idx1=RetMap(M).table_mua.Area(RetMap(M).table_mua.prf_r2>th)==1;
    idx4=RetMap(M).table_mua.Area(RetMap(M).table_mua.prf_r2>th)==4;
    
    
    %remove weird array
    idx_rm=(ecc>7 & sd<1.0);
    ecc(idx_rm)=[];
    sd(idx_rm)=[];
    idx1(idx_rm)=[];
    idx4(idx_rm)=[];
    
    v1_eccsz = [v1_eccsz; ecc(idx1) sd(idx1)];
    v4_eccsz = [v4_eccsz; ecc(idx4) sd(idx4)];
    scatter(ecc(idx1),sd(idx1));
    scatter(ecc(idx4),sd(idx4));
    
    x=[0 25];
    p1=polyfit(ecc(idx1),sd(idx1),1);
    y=polyval(p1,x);
    plot(x,y);
    p4=polyfit(ecc(idx4),sd(idx4),1);
    y=polyval(p4,x);
    plot(x,y);
    
end

%%
figure;
%subplot(1,2,1);hold on
%scatter(v1_eccsz(:,1),v1_eccsz(:,2));
%scatter(v4_eccsz(:,1),v4_eccsz(:,2));

%subplot(1,2,2);
hold on
x=[0 25];
p1=polyfit(v1_eccsz(:,1),v1_eccsz(:,2),1);
y=polyval(p1,x);
plot(x,y);
p4=polyfit(v4_eccsz(:,1),v4_eccsz(:,2),1);
y=polyval(p4,x);
plot(x,y);

poly_MUA_V1=p1;
poly_MUA_V4=p4;

% binned averages
msz1=[];
msz4=[];
step=.5:1:15.5;
for s=step
    msz1=[msz1; median(v1_eccsz(v1_eccsz(:,1)>s-0.5 & v1_eccsz(:,1)<s+0.5,2))];
    msz4=[msz4; median(v4_eccsz(v4_eccsz(:,1)>s-0.5 & v4_eccsz(:,1)<s+0.5,2))];
end

plot(step,msz1,'o')
plot(step,msz4,'o')

%%
figure;
hold on
plot(step,msz1,'o')
plot(step,msz4,'o')
x=[0 25];
p1=polyfit(step(~isnan(msz1))',msz1(~isnan(msz1)),1);
y=polyval(p1,x);
plot(x,y);
p4=polyfit(step(~isnan(msz4))',msz4(~isnan(msz4)),1);
y=polyval(p4,x);
plot(x,y);

%% example traces max neg alpha fit
M=1;
% max r2 with negative beta for alpha
maxnegbeta=max(RetMap(M).table_lfp.alpha_prf_r2(RetMap(M).table_lfp.alpha_prf_beta<0));
idx=find(RetMap(M).table_lfp.alpha_prf_r2==maxnegbeta)


% mua
muapred = RetMap(M).MUA(RetMap(M).table_mua.Instance(idx)).pred(:,RetMap(M).table_mua.ChanIndex(idx));
mua = RetMap(M).MUA(RetMap(M).table_mua.Instance(idx)).MUAm(:,RetMap(M).table_mua.ChanIndex(idx));
% lfp (same channel)
for f=1:5
    lfppred{f} = RetMap(M).LFP(RetMap(M).table_mua.Instance(idx),f).pred(:,RetMap(M).table_mua.ChanIndex(idx));
    lfp{f} = RetMap(M).LFP(RetMap(M).table_mua.Instance(idx),f).LFPm(:,RetMap(M).table_mua.ChanIndex(idx));
end

figure; 
subplot(6,1,1);hold on;
plot(muapred);
plot(mua, 'o');

for f=1:5
    subplot(6,1,1+f); hold on;
    plot(lfppred{f});
    plot(lfp{f}, 'o');
end

%% example traces max neg alpha fit
M=1;
% max r2 with positive beta for alpha
maxposbeta=max(RetMap(M).table_lfp.alpha_prf_r2(RetMap(M).table_lfp.alpha_prf_beta>0));
idx=find(RetMap(M).table_lfp.alpha_prf_r2==maxposbeta);


% mua
muapred = RetMap(M).MUA(RetMap(M).table_mua.Instance(idx)).pred(:,RetMap(M).table_mua.ChanIndex(idx));
mua = RetMap(M).MUA(RetMap(M).table_mua.Instance(idx)).MUAm(:,RetMap(M).table_mua.ChanIndex(idx));
% lfp (same channel)
for f=1:5
    lfppred{f} = RetMap(M).LFP(RetMap(M).table_mua.Instance(idx),f).pred(:,RetMap(M).table_mua.ChanIndex(idx));
    lfp{f} = RetMap(M).LFP(RetMap(M).table_mua.Instance(idx),f).LFPm(:,RetMap(M).table_mua.ChanIndex(idx));
end

figure; 
subplot(6,1,1);hold on;
plot(muapred);
plot(mua, 'o');

for f=1:5
    subplot(6,1,1+f); hold on;
    plot(lfppred{f});
    plot(lfp{f}, 'o');
end

%% are there any good negative beta-fits?
for M=1:2
    th=0.7;
    n_negtheta=sum(RetMap(M).table_lfp.theta_prf_beta<0 & RetMap(M).table_lfp.theta_prf_r2>th)
    n_negalpha=sum(RetMap(M).table_lfp.alpha_prf_beta<0 & RetMap(M).table_lfp.alpha_prf_r2>th)
    n_negbeta=sum(RetMap(M).table_lfp.beta_prf_beta<0 & RetMap(M).table_lfp.beta_prf_r2>th)
    n_lowgamma=sum(RetMap(M).table_lfp.lowgamma_prf_beta<0 & RetMap(M).table_lfp.lowgamma_prf_r2>th)
    n_highgamma=sum(RetMap(M).table_lfp.highgamma_prf_beta<0 & RetMap(M).table_lfp.highgamma_prf_r2>th)
end

%%
th=0.7;
table_lfp=[RetMap(1).table_lfp;RetMap(2).table_lfp];
table_mua=[RetMap(1).table_mua;RetMap(2).table_mua];

n_negtheta=sum(table_lfp.theta_prf_beta<0 & table_lfp.theta_prf_r2>th)
sigtheta = sum(table_lfp.theta_prf_r2>th)

n_negalpha=sum(table_lfp.alpha_prf_beta<0 & table_lfp.alpha_prf_r2>th)
sigalpha = sum(table_lfp.alpha_prf_r2>th)

n_negbeta=sum(table_lfp.beta_prf_beta<0 & table_lfp.beta_prf_r2>th)
sigbeta = sum(table_lfp.beta_prf_r2>th)

n_lowgamma=sum(table_lfp.lowgamma_prf_beta<0 & table_lfp.lowgamma_prf_r2>th)
siglowgamma = sum(table_lfp.lowgamma_prf_r2>th)

n_highgamma=sum(table_lfp.highgamma_prf_beta<0 & table_lfp.highgamma_prf_r2>th)
sighighgamma = sum(table_lfp.highgamma_prf_r2>th)

occ=[n_negtheta sigtheta;...
    n_negalpha sigalpha;...
    n_negbeta sigbeta;...
    n_lowgamma siglowgamma; ...
    n_highgamma sighighgamma];
occ_perc=(occ(:,1)./occ(:,2)).*100

sigalpha = sum(RetMap(1).table_lfp.alpha_prf_r2>th)
sigalpha = sum(RetMap(2).table_lfp.alpha_prf_r2>th)

%(table_lfp.alpha_prf_beta<0 & table_lfp.alpha_prf_r2>th)

negidx=table_lfp.alpha_prf_beta<0;
posidx=table_lfp.alpha_prf_beta>0;

rneg=table_lfp.alpha_prf_r2(negidx);
rpos=table_lfp.alpha_prf_r2(posidx);

sdneg=table_lfp.alpha_prf_sd(negidx);
sdpos=table_lfp.alpha_prf_sd(posidx);

sdneg_mua=table_mua.prf_sd(negidx);
sdneg_hg=table_lfp.highgamma_prf_sd(negidx);
sdpos_mua=table_mua.prf_sd(posidx);
sdpos_hg=table_lfp.highgamma_prf_sd(posidx);

eccneg=table_lfp.alpha_prf_ecc(negidx);
eccpos=table_lfp.alpha_prf_ecc(posidx);

areaneg=table_lfp.Area(negidx & table_lfp.alpha_prf_r2>th);
areapos=table_lfp.Area(posidx & table_lfp.alpha_prf_r2>th);

figure;
subplot(3,2,1);hold on;
hist(rpos,0:0.05:1)
title('R2 posbeta')
set(gca, 'Xlim',[0 1])
subplot(3,2,2);hold on;
hist(rneg,0:0.05:1)
title('R2 negbeta')
set(gca, 'Xlim',[0 1])

th=0.6;
npos=sum(rpos>th);
nneg=sum(rneg>th);

mSz_pos=mean(sdpos(rpos>th));
mSz_neg=mean(sdneg(rneg>th));

mEcc_pos=mean(eccpos(rpos>th));
mEcc_neg=mean(eccneg(rneg>th));

mArea_pos=mean(areapos);
mArea_neg=mean(areaneg);

subplot(3,2,3);hold on;
hist(sdpos(rpos>th),0:0.5:8)
title('SD posbeta')
set(gca, 'Xlim',[0 8.5])

subplot(3,2,4);hold on;
hist(sdneg(rneg>th),0:0.5:8)
title('SD negbeta')
set(gca, 'Xlim',[0 8.5])

subplot(3,2,5);hold on;
hist(eccpos(rpos>th),0:0.5:12)
title('ECC negbeta')
set(gca, 'Xlim',[0 12])

subplot(3,2,6);hold on;
hist(eccneg(rneg>th),0:0.5:12)
title('ECC negbeta')
set(gca, 'Xlim',[0 12])

%% plot the size of pos/neg alpha prfs against mua/high-gamma
figure; 
subplot(2,2,1);hold on
plot([0 5],[0 5],'k','LineWidth',2);
plot(sdneg_mua(rneg>0.7),sdneg(rneg>0.7),'o')
subplot(2,2,2);hold on
plot([0 5],[0 5],'k','LineWidth',2);
plot(sdneg_hg(rneg>0.7),sdneg(rneg>0.7),'o')
subplot(2,2,3);hold on
plot([0 5],[0 5],'k','LineWidth',2);
plot(sdpos_mua(rpos>0.7),sdpos(rpos>0.7),'o')
subplot(2,2,4);hold on
plot([0 5],[0 5],'k','LineWidth',2);
plot(sdpos_hg(rpos>0.7),sdpos(rpos>0.7),'o')

% t-test on pRF size comparison high gamma and pos/neg alpha
[h,p]=ttest(sdneg_hg(rneg>0.7),sdneg(rneg>0.7))
[h,p]=ttest(sdpos_hg(rpos>0.7),sdpos(rpos>0.7))

%% ECC vs Size for the LFP
th=0.7;
figure;

for a=[1 4]
    %sum(table_lfp.alpha_prf_r2>th)
    %sum(table_lfp.beta_prf_r2>th)
    
    subplot(2,2,1); hold on;
    idx1 = table_lfp.alpha_prf_r2>th & table_lfp.Area==a;
    ecc=table_lfp.alpha_prf_ecc(idx1);
    sd=table_lfp.alpha_prf_sd(idx1);
    
    if a==1
        ecc(sd>3.5)=[];
        sd(sd>3.5)=[];
    end
    
    msz1=[];
    step=.5:1:19.5;
    for s=step
        msz1=[msz1; median(sd(ecc>s-0.5 & ecc<s+0.5))];
    end

    x=[0 20];
    %p1=polyfit(ecc,sd,1);
    %y=polyval(p1,x);
    %plot(ecc,sd,'o',x,y);
    
    p1=polyfit(step(~isnan(msz1)),msz1(~isnan(msz1))',1);
    y=polyval(p1,x);
    plot(step,msz1,'o',x,y);
    
    if a==1
        poly_LFPa_V1 = p1;
    elseif a==4
        poly_LFPa_V4 = p1;
    end
    
    subplot(2,2,2); hold on;
    idx1 = table_lfp.beta_prf_r2>th & table_lfp.Area==a;
    ecc=table_lfp.beta_prf_ecc(idx1);
    sd=table_lfp.beta_prf_sd(idx1);
    
    if a==1
        ecc(sd>3.5)=[];
        sd(sd>3.5)=[];
    end
    
    msz1=[];
    for s=step
        msz1=[msz1; median(sd(ecc>s-0.5 & ecc<s+0.5))];
    end

    x=[0 20];
    %p1=polyfit(ecc,sd,1);
    %y=polyval(p1,x);
    %plot(ecc,sd,'o',x,y);
    
    p1=polyfit(step(~isnan(msz1)),msz1(~isnan(msz1))',1);
    y=polyval(p1,x);
    plot(step,msz1,'o',x,y);
    
    if a==1
        poly_LFPb_V1 = p1;
    elseif a==4
        poly_LFPb_V4 = p1;
    end
    
    subplot(2,2,3); hold on;
    idx1 = table_lfp.lowgamma_prf_r2>th & table_lfp.Area==a;
    ecc=table_lfp.lowgamma_prf_ecc(idx1);
    sd=table_lfp.lowgamma_prf_sd(idx1);
    
    if a==1
        ecc(sd>3.5)=[];
        sd(sd>3.5)=[];
    end
    
    msz1=[];
    for s=step
        msz1=[msz1; median(sd(ecc>s-0.5 & ecc<s+0.5))];
    end

    x=[0 20];
    %p1=polyfit(ecc,sd,1);
    %y=polyval(p1,x);
    %plot(ecc,sd,'o',x,y);
    
    p1=polyfit(step(~isnan(msz1)),msz1(~isnan(msz1))',1);
    y=polyval(p1,x);
    plot(step,msz1,'o',x,y);
    
    if a==1
        poly_LFPg1_V1 = p1;
    elseif a==4
        poly_LFPg1_V4 = p1;
    end
    
    subplot(2,2,4); hold on;
    idx1 = table_lfp.highgamma_prf_r2>th & table_lfp.Area==a;
    ecc=table_lfp.highgamma_prf_ecc(idx1);
    sd=table_lfp.highgamma_prf_sd(idx1);
    
    if a==1
        ecc(sd>3.5)=[];
        sd(sd>3.5)=[];
    end
    
    msz1=[];
    for s=step
        msz1=[msz1; median(sd(ecc>s-0.5 & ecc<s+0.5))];
    end

    x=[0 20];
    %p1=polyfit(ecc,sd,1);
    %y=polyval(p1,x);
    %plot(ecc,sd,'o',x,y);
    
    p1=polyfit(step(~isnan(msz1)),msz1(~isnan(msz1))',1);
    y=polyval(p1,x);
    plot(step,msz1,'o',x,y);
    
    if a==1
        poly_LFPg2_V1 = p1;
    elseif a==4
        poly_LFPg2_V4 = p1;
    end
end

%% crossmodal comparisons
th=0.7;
table_mua=[RetMap(1).table_mua;RetMap(2).table_mua];

figure;
cc=zeros(4,4);pp=zeros(4,4);
% R2 based on MUA

subplot(4,4,1);hold on;
idx1=table_mua.prf_r2>th;
idx2=table_lfp.alpha_prf_r2>th;
idx3=logical(idx1.*idx2);

plot([0 8],[0 8]); hold on;
plot(table_mua.prf_sd(idx3),table_lfp.alpha_prf_sd(idx3),'o')
[r,p]=corrcoef(table_mua.prf_sd(idx3),table_lfp.alpha_prf_sd(idx3));
cc(1,1)=r(2);pp(1,1)=p(2);

nAlpha_sz = table_lfp.alpha_prf_sd(idx3)./table_mua.prf_sd(idx3);

[x1,y1]=pol2cart(...
    deg2rad(table_mua.prf_pol(idx3)), deg2rad(table_mua.prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.alpha_prf_pol(idx3)), deg2rad(table_lfp.alpha_prf_ecc(idx3)) ...
    );
dAlphaMUA=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,5);hold on;
idx1=table_mua.prf_r2>th;
idx2=table_lfp.beta_prf_r2>th;
idx2b=table_lfp.beta_prf_sd<8;
idx3=logical(idx1.*idx2.*idx2b);

plot([0 8],[0 8]); hold on;
plot(table_mua.prf_sd(idx3),table_lfp.beta_prf_sd(idx3),'o')
[r,p]=corrcoef(table_mua.prf_sd(idx3),table_lfp.beta_prf_sd(idx3));
cc(2,1)=r(2);pp(2,1)=p(2);

nBeta_sz = table_lfp.beta_prf_sd(idx3)./table_mua.prf_sd(idx3);

[x1,y1]=pol2cart(...
    deg2rad(table_mua.prf_pol(idx3)), deg2rad(table_mua.prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.beta_prf_pol(idx3)), deg2rad(table_lfp.beta_prf_ecc(idx3)) ...
    );
dBetaMUA=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,9);hold on;
idx1=table_mua.prf_r2>th;
idx2=table_lfp.lowgamma_prf_r2>th;
idx3=logical(idx1.*idx2);

plot([0 8],[0 8]); hold on;
plot(table_mua.prf_sd(idx3),table_lfp.lowgamma_prf_sd(idx3),'o')
[r,p]=corrcoef(table_mua.prf_sd(idx3),table_lfp.lowgamma_prf_sd(idx3));
cc(3,1)=r(2);pp(3,1)=p(2);

nGamma1_sz = table_lfp.lowgamma_prf_sd(idx3)./table_mua.prf_sd(idx3);

[x1,y1]=pol2cart(...
    deg2rad(table_mua.prf_pol(idx3)), deg2rad(table_mua.prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.lowgamma_prf_pol(idx3)), deg2rad(table_lfp.lowgamma_prf_ecc(idx3)) ...
    );
dGam1MUA=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,13);hold on;
idx1=table_mua.prf_r2>th;
idx2=table_lfp.highgamma_prf_r2>th;
idx3=logical(idx1.*idx2);

plot([0 8],[0 8]); hold on;
plot(table_mua.prf_sd(idx3),table_lfp.highgamma_prf_sd(idx3),'o')
[r,p]=corrcoef(table_mua.prf_sd(idx3),table_lfp.highgamma_prf_sd(idx3));
cc(4,1)=r(2);pp(4,1)=p(2);

nGamma2_sz = table_lfp.highgamma_prf_sd(idx3)./table_mua.prf_sd(idx3);

[x1,y1]=pol2cart(...
    deg2rad(table_mua.prf_pol(idx3)), deg2rad(table_mua.prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.highgamma_prf_pol(idx3)), deg2rad(table_lfp.highgamma_prf_ecc(idx3)) ...
    );
dGam2MUA=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,6);hold on;
idx1=table_lfp.alpha_prf_r2>th;
idx2=table_lfp.beta_prf_r2>th;
idx3=logical(idx1.*idx2);

plot([0 8],[0 8]); hold on;
plot(table_lfp.alpha_prf_sd(idx3),table_lfp.beta_prf_sd(idx3),'o')
[r,p]=corrcoef(table_lfp.alpha_prf_sd(idx3),table_lfp.beta_prf_sd(idx3));
cc(2,2)=r(2); pp(2,2)=p(2);

[x1,y1]=pol2cart(...
    deg2rad(table_lfp.alpha_prf_pol(idx3)), deg2rad(table_lfp.alpha_prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.beta_prf_pol(idx3)), deg2rad(table_lfp.beta_prf_ecc(idx3)) ...
    );
dAlphaBeta=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,10);hold on;
idx1=table_lfp.alpha_prf_r2>th;
idx2=table_lfp.lowgamma_prf_r2>th;
idx2b=table_lfp.lowgamma_prf_sd<8;
idx3=logical(idx1.*idx2.*idx2b);

plot([0 8],[0 8]); hold on;
plot(table_lfp.alpha_prf_sd(idx3),table_lfp.lowgamma_prf_sd(idx3),'o')
[r,p]=corrcoef(table_lfp.alpha_prf_sd(idx3),table_lfp.lowgamma_prf_sd(idx3));
cc(3,2)=r(2);pp(3,2)=p(2);

[x1,y1]=pol2cart(...
    deg2rad(table_lfp.alpha_prf_pol(idx3)), deg2rad(table_lfp.alpha_prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.lowgamma_prf_pol(idx3)), deg2rad(table_lfp.lowgamma_prf_ecc(idx3)) ...
    );
dAlphaGam1=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,14);hold on;
idx1=table_lfp.alpha_prf_r2>th;
idx2=table_lfp.highgamma_prf_r2>th;
idx2b=table_lfp.highgamma_prf_sd<8;
idx3=logical(idx1.*idx2.*idx2b);

plot([0 8],[0 8]); hold on;
plot(table_lfp.alpha_prf_sd(idx3),table_lfp.lowgamma_prf_sd(idx3),'o')
[r,p]=corrcoef(table_lfp.alpha_prf_sd(idx3),table_lfp.highgamma_prf_sd(idx3));
cc(4,2)=r(2);pp(4,1)=p(2);

[x1,y1]=pol2cart(...
    deg2rad(table_lfp.alpha_prf_pol(idx3)), deg2rad(table_lfp.alpha_prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.highgamma_prf_pol(idx3)), deg2rad(table_lfp.highgamma_prf_ecc(idx3)) ...
    );
dAlphaGam2=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,11);hold on;
idx1=table_lfp.beta_prf_r2>th;
idx2=table_lfp.lowgamma_prf_r2>th;
idx3=logical(idx1.*idx2);

plot([0 8],[0 8]); hold on;
plot(table_lfp.beta_prf_sd(idx3),table_lfp.lowgamma_prf_sd(idx3),'o')
[r,p]=corrcoef(table_lfp.beta_prf_sd(idx3),table_lfp.lowgamma_prf_sd(idx3));
cc(3,3)=r(2);pp(3,3)=p(2);

[x1,y1]=pol2cart(...
    deg2rad(table_lfp.beta_prf_pol(idx3)), deg2rad(table_lfp.beta_prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.lowgamma_prf_pol(idx3)), deg2rad(table_lfp.lowgamma_prf_ecc(idx3)) ...
    );
dBetaGam1=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,15);hold on;
idx1=table_lfp.beta_prf_r2>th;
idx2=table_lfp.highgamma_prf_r2>th;
idx2b=table_lfp.highgamma_prf_sd<8;
idx3=logical(idx1.*idx2.*idx2b);
plot([0 8],[0 8]); hold on;
plot(table_lfp.beta_prf_sd(idx3),table_lfp.highgamma_prf_sd(idx3),'o')
[r,p]=corrcoef(table_lfp.beta_prf_sd(idx3),table_lfp.highgamma_prf_sd(idx3));
cc(3,4)=r(2);pp(3,4)=p(2);

[x1,y1]=pol2cart(...
    deg2rad(table_lfp.beta_prf_pol(idx3)), deg2rad(table_lfp.beta_prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.highgamma_prf_pol(idx3)), deg2rad(table_lfp.highgamma_prf_ecc(idx3)) ...
    );
dBetaGam2=sqrt(((x2-x1).^2)+((y2-y1).^2));

subplot(4,4,16);hold on;
idx1=table_lfp.lowgamma_prf_r2>th;
idx2=table_lfp.highgamma_prf_r2>th;
idx3=logical(idx1.*idx2.*idx2b);

plot([0 8],[0 8]); hold on;
plot(table_lfp.lowgamma_prf_sd(idx3),table_lfp.highgamma_prf_sd(idx3),'o')
[r,p]=corrcoef(table_lfp.lowgamma_prf_sd(idx3),table_lfp.highgamma_prf_sd(idx3));
cc(4,4)=r(2);pp(4,4)=p(2);

[x1,y1]=pol2cart(...
    deg2rad(table_lfp.lowgamma_prf_pol(idx3)), deg2rad(table_lfp.lowgamma_prf_ecc(idx3)) ...
    );
[x2,y2]=pol2cart(...
    deg2rad(table_lfp.highgamma_prf_pol(idx3)), deg2rad(table_lfp.highgamma_prf_ecc(idx3)) ...
    );
dGam1Gam2=sqrt(((x2-x1).^2)+((y2-y1).^2));

mFreq = [mean(nAlpha_sz) mean(nBeta_sz) mean(nGamma1_sz) mean(nGamma2_sz)];
seFreq = [std(nAlpha_sz)./sqrt(length(nAlpha_sz)) ...
    std(nBeta_sz)./sqrt(length(nBeta_sz)) ...
    std(nGamma1_sz)./sqrt(length(nGamma1_sz)) ...
    std(nGamma2_sz)./sqrt(length(nGamma2_sz))];
 
 mD=[mean(dAlphaMUA) 0 0 0;...
     mean(dBetaMUA) mean(dAlphaBeta) 0 0;...
     mean(dGam1MUA) mean(dAlphaGam1) mean(dBetaGam1) 0;...
     mean(dAlphaMUA) mean(dAlphaGam2) mean(dBetaGam2) mean(dGam1Gam2)];
 
 %%

