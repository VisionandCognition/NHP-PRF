%% stats for linear models comparison
% do V1 only for now

%% DATA
ecc=[]; sz=[]; % eccentricity and size  
signal = []; % categorical 1=MRI, 2=MUA, 3-7=LFP bands
area=[]; % categorical 1=V1,4=V4
model=[]; % categorical 1=lin, 2=lin_ng, 3=css, 4=dog

% v1
ecc = [ecc; ...
    mri1(1).ECC; ...
    mua1(1).ECC; ...
     ];
sz = [ sz; ...
    mri1(1).S; ...
    mua1(1).S; ...
    ];
signal = [signal;...
    1*ones(size(mri1(1).S)); ...
    2*ones(size(mua1(1).S)); ...
    ];
area = [area; ...
    1*ones(size(mri1(1).S)); ...
    1*ones(size(mua1(1).S)); ...
    ];


for i=1:5
    ecc = [ecc; lfp1(i).ECC];
    sz = [ sz; lfp1(i).S];
    signal = [signal; (2+i)*ones(size(lfp1(i).S))];
    area = [area; 1*ones(size(lfp1(i).S))];
end

% v4
ecc = [ecc; ...
    mri4(1).ECC; ...
    mua4(1).ECC; ...
     ];
sz = [ sz; ...
    mri4(1).S; ...
    mua4(1).S; ...
    ];
signal = [signal;...
    1*ones(size(mri4(1).S)); ...
    2*ones(size(mua4(1).S)); ...
    ];
area = [area; ...
    4*ones(size(mri4(1).S)); ...
    4*ones(size(mua4(1).S)); ...
    ];


for i=1:5
    ecc = [ecc; lfp4(i).ECC];
    sz = [ sz; lfp4(i).S];
    signal = [signal; (2+i)*ones(size(lfp4(i).S))];
    area = [area; 4*ones(size(lfp4(i).S))];
end

%
sigtype = {1,2,3,4,5,6,7;'mri' ,'mua','alpha','beta','theta','hgam','lgam'};
signal=categorical(signal);
area=categorical(area);
areaname = {' V1', 'V4'};
DT=table(signal,area, ecc, sz);

%%
figure;

% do the modelfits
DT1=DT(DT.area==categorical(1) & DT.ecc<=MaxECC,:);
mdl_v1 = fitlm(DT1,'sz ~ 1 + signal*ecc');
subplot(1,2,1); gscatter(DT1.ecc,DT1.sz,DT1.signal,cc(1:3,:),'.',10)

DT4=DT(DT.area==categorical(4) & DT.ecc<=MaxECC,:);
mdl_v4 = fitlm(DT4,'sz ~ 1 + signal*ecc');
subplot(1,2,2); gscatter(DT4.ecc,DT4.sz,DT4.signal,cc(1:3,:),'.',10)

% Look at effect sizes ===
%plotEffects(mdl_v1)
%plotInteraction(mdl_v1,'ecc','signal')

%% V1 ---
display(mdl_v1);
display(mdl_v1.CoefficientNames);
display(sigtype);

% stats structure
stats.V1.mdl = mdl_v1;
stats.V1.signals = sigtype;

SigToCheck = [];
for i=2:length(mdl_v1.CoefficientNames)/2
    SigToCheck = [SigToCheck ...
        str2num(mdl_v1.CoefficientNames{i}(8:end))];
end

stats.V1.sig(1).name = 'mri';
stats.V1.sig(1).ic.val = mdl_v1.Coefficients.Estimate(1);
stats.V1.sig(1).ic.SE = mdl_v1.Coefficients.SE(1);
stats.V1.sig(1).ic.T = mdl_v1.Coefficients.tStat(1);
stats.V1.sig(1).ic.p = mdl_v1.Coefficients.pValue(1);

stats.V1.sig(1).sl.val = mdl_v1.Coefficients.Estimate(1+length(mdl_v1.CoefficientNames)/2);
stats.V1.sig(1).sl.SE = mdl_v1.Coefficients.SE(1+length(mdl_v1.CoefficientNames)/2);
stats.V1.sig(1).sl.T = mdl_v1.Coefficients.tStat(1+length(mdl_v1.CoefficientNames)/2);
stats.V1.sig(1).sl.p = mdl_v1.Coefficients.pValue(1+length(mdl_v1.CoefficientNames)/2);

for ST = 1:length(SigToCheck)
    stats.V1.sig(1+ST).name = sigtype{2,SigToCheck(ST)};
    stats.V1.sig(1+ST).ic.val = stats.V1.sig(1).ic.val + ...
        mdl_v1.Coefficients.Estimate(1+ST);
    stats.V1.sig(1+ST).ic.SE = mdl_v1.Coefficients.SE(1+ST);
    stats.V1.sig(1+ST).ic.T = mdl_v1.Coefficients.tStat(1+ST);
    stats.V1.sig(1+ST).ic.p = mdl_v1.Coefficients.pValue(1+ST);

    stats.V1.sig(1+ST).sl.val = stats.V1.sig(1).sl.val + ...
        mdl_v1.Coefficients.Estimate(ST+1+length(mdl_v1.CoefficientNames)/2);
    stats.V1.sig(1+ST).sl.SE = mdl_v1.Coefficients.SE(ST+1+length(mdl_v1.CoefficientNames)/2);
    stats.V1.sig(1+ST).sl.T = mdl_v1.Coefficients.tStat(ST+1+length(mdl_v1.CoefficientNames)/2);
    stats.V1.sig(1+ST).sl.p = mdl_v1.Coefficients.pValue(ST+1+length(mdl_v1.CoefficientNames)/2);
    
    % intercept
    cm=zeros(1,length(mdl_v1.CoefficientNames)); 
    cm(1)=-1;
    cm(1+ST)=1;
        
    [p,f,r] = coefTest(mdl_v1,cm);
    stats.V1.SigVsMRI(ST).sig = sigtype{2,SigToCheck(ST)};
    stats.V1.SigVsMRI(ST).ic.p = p; 
    stats.V1.SigVsMRI(ST).ic.f = f; 
    stats.V1.SigVsMRI(ST).ic.r = r;
    
    % slope
    cm=zeros(1,length(mdl_v1.CoefficientNames)); 
    cm(1+length(mdl_v1.CoefficientNames)/2)=-1;
    cm(1+(length(mdl_v1.CoefficientNames)/2)+ST)=1;
    [p,f,r] = coefTest(mdl_v1,cm);
    stats.V1.SigVsMRI(ST).sl.p = p; 
    stats.V1.SigVsMRI(ST).sl.f = f; 
    stats.V1.SigVsMRI(ST).sl.r = r;

    % intercept
    cm=zeros(1,length(mdl_v1.CoefficientNames)); 
    cm([1,1+length(mdl_v1.CoefficientNames)/2])=-1;
    cm([1+ST,(1+length(mdl_v1.CoefficientNames)/2)+ST])=1;
    [p,f,r] = coefTest(mdl_v1,cm);
    stats.V1.SigVsMRI(ST).icsl.p = p; 
    stats.V1.SigVsMRI(ST).icsl.f = f; 
    stats.V1.SigVsMRI(ST).icsl.r = r;
end

% collect
stats.V1.SigVsMRI_all={'sig','icf','icp','slf','slp'};
for i=1:length(stats.V1.SigVsMRI)
    stats.V1.SigVsMRI_all = [...
        stats.V1.SigVsMRI_all; {...
        stats.V1.SigVsMRI(i).sig, ...
        stats.V1.SigVsMRI(i).ic.f, ...
        stats.V1.SigVsMRI(i).ic.p,...
        stats.V1.SigVsMRI(i).sl.f, ...
        stats.V1.SigVsMRI(i).sl.p}];
end

%% V4 ---
display(mdl_v4);
display(mdl_v4.CoefficientNames);
display(sigtype);

% stats structure
stats.V4.mdl = mdl_v4;
stats.V4.signals = sigtype;

SigToCheck = [];
for i=2:length(mdl_v4.CoefficientNames)/2
    SigToCheck = [SigToCheck ...
        str2num(mdl_v4.CoefficientNames{i}(8:end))];
end

stats.V4.sig(1).name = 'mri';
stats.V4.sig(1).ic.val = mdl_v4.Coefficients.Estimate(1);
stats.V4.sig(1).ic.SE = mdl_v4.Coefficients.SE(1);
stats.V4.sig(1).ic.T = mdl_v4.Coefficients.tStat(1);
stats.V4.sig(1).ic.p = mdl_v4.Coefficients.pValue(1);

stats.V4.sig(1).sl.val = mdl_v4.Coefficients.Estimate(1+length(mdl_v4.CoefficientNames)/2);
stats.V4.sig(1).sl.SE = mdl_v4.Coefficients.SE(1+length(mdl_v4.CoefficientNames)/2);
stats.V4.sig(1).sl.T = mdl_v4.Coefficients.tStat(1+length(mdl_v4.CoefficientNames)/2);
stats.V4.sig(1).sl.p = mdl_v4.Coefficients.pValue(1+length(mdl_v4.CoefficientNames)/2);

for ST = 1:length(SigToCheck)
    stats.V4.sig(1+ST).name = sigtype{2,SigToCheck(ST)};
    stats.V4.sig(1+ST).ic.val = stats.V4.sig(1).ic.val + ...
        mdl_v4.Coefficients.Estimate(1+ST);
    stats.V4.sig(1+ST).ic.SE = mdl_v4.Coefficients.SE(1+ST);
    stats.V4.sig(1+ST).ic.T = mdl_v4.Coefficients.tStat(1+ST);
    stats.V4.sig(1+ST).ic.p = mdl_v4.Coefficients.pValue(1+ST);

    stats.V4.sig(1+ST).sl.val = stats.V4.sig(1).sl.val + ...
        mdl_v4.Coefficients.Estimate(ST+1+length(mdl_v4.CoefficientNames)/2);
    stats.V4.sig(1+ST).sl.SE = mdl_v4.Coefficients.SE(ST+1+length(mdl_v4.CoefficientNames)/2);
    stats.V4.sig(1+ST).sl.T = mdl_v4.Coefficients.tStat(ST+1+length(mdl_v4.CoefficientNames)/2);
    stats.V4.sig(1+ST).sl.p = mdl_v4.Coefficients.pValue(ST+1+length(mdl_v4.CoefficientNames)/2);
    
    % intercept
    cm=zeros(1,length(mdl_v4.CoefficientNames)); 
    cm(1)=-1;
    cm(1+ST)=1;
        
    [p,f,r] = coefTest(mdl_v4,cm);
    stats.V4.SigVsMRI(ST).sig = sigtype{2,SigToCheck(ST)};
    stats.V4.SigVsMRI(ST).ic.p = p; 
    stats.V4.SigVsMRI(ST).ic.f = f; 
    stats.V4.SigVsMRI(ST).ic.r = r;
    
    % slope
    cm=zeros(1,length(mdl_v4.CoefficientNames)); 
    cm(1+length(mdl_v4.CoefficientNames)/2)=-1;
    cm(1+(length(mdl_v4.CoefficientNames)/2)+ST)=1;
    [p,f,r] = coefTest(mdl_v4,cm);
    stats.V4.SigVsMRI(ST).sl.p = p; 
    stats.V4.SigVsMRI(ST).sl.f = f; 
    stats.V4.SigVsMRI(ST).sl.r = r;

    % intercept
    cm=zeros(1,length(mdl_v4.CoefficientNames)); 
    cm([1,1+length(mdl_v4.CoefficientNames)/2])=-1;
    cm([1+ST,(1+length(mdl_v4.CoefficientNames)/2)+ST])=1;
    [p,f,r] = coefTest(mdl_v4,cm);
    stats.V4.SigVsMRI(ST).icsl.p = p; 
    stats.V4.SigVsMRI(ST).icsl.f = f; 
    stats.V4.SigVsMRI(ST).icsl.r = r;
end





%% fitlme (with 'monkey')??