%% Set-up =================================================================
% figure saving folder
pngfld = fullfile(pwd,'fig_png');
svgfld = fullfile(pwd,'fig_svg');
[~,~] = mkdir(pngfld);
[~,~] = mkdir(svgfld);

% load result if it exists
if isfile('ExampleElectrode.mat')
    LoadedResults=true;
    load('ExampleElectrode.mat');
else
    LoadedResults=false;
end

addpath(genpath('/Users/chris/Documents/MRI_ANALYSIS/NHP-PRF/prfCode/LISA/OnServer/analyzePRF'))
modeltype = 'linear_ephys'; TR=0.5;

%% Load files =============================================================
if ~LoadedResults
    RESULT_FILE = [ '/home/chris/Documents/MRI_ANALYSIS/NHP-PRF/FitResults/'...
        'ephys/Combined/ORG/AllFits_ephys_cv1.mat'];
    load(RESULT_FILE);
    RES = R(1).model(1);
    
    maxi=[];
    for i=1:8
        maxi = [maxi; i max(RES.MUA(i).R2(:,1)) ...
            find(RES.MUA(i).R2(:,1)==max(RES.MUA(i).R2(:,1)),1,'first')];
    end
    lidx = find(maxi(:,2)==max(maxi(:,2)),1,'first');
    inst_elec = [lidx max(lidx,3)];
    inst=lidx; ch=inst_elec(2);
    
    load(['/home/chris/Documents/MRI_ANALYSIS/NHP-PRF/Data/ephys/lick/mua/'...
        'Lick_20180807_B2_array_' num2str(lidx) '_mMUA.mat']);
    %mMUA and mLFP hold the data
    load(['/home/chris/Documents/MRI_ANALYSIS/NHP-PRF/Data/ephys/lick/lfp/'...
        'Lick_20180807_B2_array_' num2str(lidx) '_mLFP.mat']);
    warning off;
    load('/home/chris/Documents/MRI_ANALYSIS/NHP-PRF/Data/ephys/lick/stim/Lick_STIM.mat');
    warning on;
    cv=0;
end

%% Prep data and run modelfit =============================================
if ~LoadedResults
    mua_data={};
    fprintf('Concatenating stimuli and volumes...\n');
    mua_data{1} = []; mua_data{1} = cat(1,mua_data{1},mMUA(ch).bar);
    stimulus{1}=[];
    for imgnr=1:length(STIM.img)
        stimulus{1}=cat(3,stimulus{1},STIM.img{imgnr});
    end
    
    inv_idx = [150:-1:121 180:-1:151 210:-1:181 240:-1:211];
    mua_data2{1}=[];
    
    for elec = 1:size(mua_data{1},1)
        mua_data2{1} = cat(1,mua_data2{1},...
            mean([mua_data{1}(elec,1:120);mua_data{1}(elec,inv_idx)],1));
    end
    mua_data_org = mua_data; mua_data = mua_data2; clear ephys_data2
    
    stimulus{1}=[];
    for imgnr=1:length(STIM.img)
        % RESAMPLE STIMULUS >> 295 x 295 means 10px = 1 deg
        rsIMG = imresize(STIM.img{imgnr} ,[295 295]);
        stimulus{1}=cat(3,stimulus{1},rsIMG);
    end
    
    options = R(1).model(1).MUA(1).options;
    options.xvalmode = cv;
    options.maxpolydeg = 0;
    options.vxs = 1;
    
    % run modelfit
    result_mua= analyzePRF_modeltype(stimulus,mua_data_org,TR,options,modeltype);
    
    for fb=1:5 % loop over frequency bands
        % concatenate -----
        lfp_data={};
        fprintf(['Frequency band ' num2str(fb) '\n']);
        fprintf('Concatenating stimuli and volumes...\n');
        
        lfp_data{1} = [];
        lfp_data{1}=cat(1,lfp_data{1},...
            mLFP(ch).freq(fb).bar - mLFP(ch).freq(fb).BL);
        
        inv_idx = [150:-1:121 180:-1:151 210:-1:181 240:-1:211];
        lfp_data2{1}=[];
        for elec = 1:size(lfp_data{1},1)
            lfp_data2{1} = cat(1,lfp_data2{1},...
                mean([lfp_data{1}(elec,1:120);lfp_data{1}(elec,inv_idx)],1));
        end
        
        lfp_data_org{fb} = lfp_data;  lfp_data = lfp_data2; clear lfp_data2
        result_lfp(fb) = analyzePRF_modeltype(stimulus,lfp_data_org{fb},TR,options,modeltype);
    end
end

%% Save the results =======================================================
if ~LoadedResults
    save('ExampleElectrode','TR','stimulus','modeltype',...
        'mua_data_org','result_mua',...
        'lfp_data_org','result_lfp');
end

%% Plot fit prediction and data together ==================================
% === MUA ===
data = mua_data_org; tr=TR;  
res = [size(stimulus{1},1) size(stimulus{1},2)];
resmx = max(res);                   
degs = result_mua.options.maxpolydeg; 
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);
options = result_mua.options;
   
stimulusPP = {};
for p=1:length(stimulus)
    stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
    stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

switch modeltype
    case 'css_hrf'
        % -- CSS --
        % define the model (parameters are R C S G N)
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
            makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
            (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),options.hrf,dd(:,prod(res)+1));
    case 'linear_hrf'
        % -- Linear (skip second step where exponential is fit) --
        % define the model (parameters are R C S G)
        if options.allowneggain
            modelfun = @(pp,dd) conv2run(pp(4) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]),options.hrf,dd(:,prod(res)+1));
        else
            modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]),options.hrf,dd(:,prod(res)+1));

        end
    case 'dog_hrf'
        % -- DOG (center surround
        % define the model (parameters are R C S G sdratio normamp)
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res), ...
            (makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2)) - ...
            (pp(6) .* makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)/pp(5)),abs(pp(3)/pp(5)),xx,yy,0,0) / (2*pi*abs(pp(3)/pp(5))^2)) ...
            )) ; 0]),options.hrf,dd(:,prod(res)+1));
    case 'css_ephys'
        % -- CSS without convolution with HRF
        % define the model (parameters are R C S G N)
        modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
            makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
            (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5));
    case 'linear_ephys'
        % -- Linear without convolution with HRF
        % define the model (parameters are R C S G)
                
        if options.allowneggain
            modelfun = @(pp,dd) pp(4) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]);
        else
            modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]);
        end
    case 'dog_ephys'
        % -- DOG without convolution with HRF
        % define the model (parameters are R C S G sdratio normamp)
        modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res), ...
            (makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2)) - ...
            (pp(6) .* makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)/pp(5)),abs(pp(3)/pp(5)),xx,yy,0,0) / (2*pi*abs(pp(3)/pp(5))^2)) ...
            )) ; 0]);
end

polymatrix = {};
for p=1:length(degs)
    polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
end
   
% Which channel should we inspect? 
ch = 1; %4;
    
datats = {}; modelts = {};
for p=1:length(data)
    datats{p} =  polymatrix{p}*data{p}(ch,:)';
    modelts{p} = polymatrix{p}*modelfun(result_mua.params(1,:,ch),single(stimulusPP{p}));
end
    
% Visualize the results
f=figure; hold on;
set(gcf,'Units','points','Position',[100 100 1000 300]);
plot(.5:0.5:120,cat(1,datats{:}),'ok','MarkerSize',6,'MarkerFaceColor',[.75 .75 .75]);
plot(.5:0.5:120,cat(1,modelts{:}),'k','Linewidth',2);
xlabel('Time (stimulus position)'); ylabel('Signal');
ax = axis; set(gca,'xlim',[0 121]);
%axis([.5 1200+.5 ax(3:4)]);
title('Time-series data MUA - Example electrode');
legend({'Data','Model'})
saveas(f,fullfile(pngfld, 'ExampleElectrode_MUA.png'));
saveas(f,fullfile(svgfld, 'ExampleElectrode_MUA.svg'));

% === LFP ===
for fb=1:5
    data = lfp_data_org{fb};
    res = [size(stimulus{1},1) size(stimulus{1},2)];
    resmx = max(res);                   
    degs = result_lfp(fb).options.maxpolydeg; 
    [d,xx,yy] = makegaussian2d(resmx,2,2,2,2);
    options = result_lfp(fb).options;

    stimulusPP = {};
    for p=1:length(stimulus)
        stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
        stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
    end

    switch modeltype
    case 'css_hrf'
        % -- CSS --
        % define the model (parameters are R C S G N)
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
            makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
            (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),options.hrf,dd(:,prod(res)+1));
    case 'linear_hrf'
        % -- Linear (skip second step where exponential is fit) --
        % define the model (parameters are R C S G)
        if options.allowneggain
            modelfun = @(pp,dd) conv2run(pp(4) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]),options.hrf,dd(:,prod(res)+1));
        else
            modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]),options.hrf,dd(:,prod(res)+1));

        end
    case 'dog_hrf'
        % -- DOG (center surround
        % define the model (parameters are R C S G sdratio normamp)
        modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res), ...
            (makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2)) - ...
            (pp(6) .* makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)/pp(5)),abs(pp(3)/pp(5)),xx,yy,0,0) / (2*pi*abs(pp(3)/pp(5))^2)) ...
            )) ; 0]),options.hrf,dd(:,prod(res)+1));
    case 'css_ephys'
        % -- CSS without convolution with HRF
        % define the model (parameters are R C S G N)
        modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
            makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
            (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5));
    case 'linear_ephys'
        % -- Linear without convolution with HRF
        % define the model (parameters are R C S G)
                
        if options.allowneggain
            modelfun = @(pp,dd) pp(4) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]);
        else
            modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
                makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
                (2*pi*abs(pp(3))^2))); 0]);
        end
    case 'dog_ephys'
        % -- DOG without convolution with HRF
        % define the model (parameters are R C S G sdratio normamp)
        modelfun = @(pp,dd) posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res), ...
            (makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2)) - ...
            (pp(6) .* makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)/pp(5)),abs(pp(3)/pp(5)),xx,yy,0,0) / (2*pi*abs(pp(3)/pp(5))^2)) ...
            )) ; 0]);
    end

    polymatrix = {};
    for p=1:length(degs)
        polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
    end
   
    % Which channel should we inspect? 
    ch = 1; %4;
    
    datats = {}; modelts = {};
    for p=1:length(data)
        datats{p} =  polymatrix{p}*data{p}(ch,:)';
        modelts{p} = polymatrix{p}*modelfun(result_lfp(fb).params(1,:,ch),single(stimulusPP{p}));
    end
    
    % Visualize the results
    f=figure; hold on;
    set(gcf,'Units','points','Position',[100 100 1000 300]);
    plot(.5:0.5:120,cat(1,datats{:}),'ok','MarkerSize',6,'MarkerFaceColor',[.75 .75 .75]);
    plot(.5:0.5:120,cat(1,modelts{:}),'k','Linewidth',2);
    xlabel('Time (stimulus position)'); ylabel('Signal');
    ax = axis; set(gca,'xlim',[0 121]);
    %axis([.5 1200+.5 ax(3:4)]);
    title(['Time-series data LFP fb-' num2str(fb) ' - Example electrode']); 
    legend({'Data','Model'})
    saveas(f,fullfile(pngfld, ['ExampleElectrode_LFP-' num2str(fb) '.png']));
    saveas(f,fullfile(svgfld, ['ExampleElectrode_LFP-' num2str(fb) '.svg']));
end

close all;






















