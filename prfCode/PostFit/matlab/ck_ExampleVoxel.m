%% Set-up =================================================================
% figure saving folder
pngfld = fullfile(pwd,'fig_png');
svgfld = fullfile(pwd,'fig_svg');
[~,~] = mkdir(pngfld);
[~,~] = mkdir(svgfld);

% load result if it exists
if isfile('ExampleVoxel.mat')
    LoadedResults=true;
    load('ExampleVoxel.mat');
else
    LoadedResults=false;
end

addpath(genpath('/Users/chris/Documents/MRI_ANALYSIS/NHP-PRF/prfCode/LISA/OnServer/analyzePRF'))

%% Load files =============================================================
if ~LoadedResults
    RESULT_FILE = ...
        ['/home/chris/Documents/MRI_ANALYSIS/NHP-PRF/FitResults/MRI/'...
        'danny/linhrf_cv1_mhrf/pRF_Sess-linhrf_cv1_mhrf.mat'];
    load(RESULT_FILE);
    modeltype = 'linear_hrf';
    
    DATA_FILE = ...
        ['/home/chris/Documents/MRI_ANALYSIS/NHP-PRF/Data/MRI/cv/'...
        'danny/AllSessions-avg-cv.mat'];
    load(DATA_FILE,'stim');
    load(DATA_FILE,'sess_wmeanBOLD');
    load(DATA_FILE,'sess_wmeanBOLD_inv');
end

%% Prep data and run modelfit =============================================
if ~LoadedResults
    % Find voxel with highest R2
    R2_1 = result.R2(:,:,:,1);
    % get 3d coordinate idx of this voxel
    [vx,vy,vz] = ind2sub(size(R2_1), find(R2_1==max(R2_1(:)),1,'first'));
    vidx = find(R2_1==max(R2_1(:)),1,'first');
    
    % concat into shape analyzePRF wants
    for RUNNR = 1:length(sess_wmeanBOLD)
        fmri_data{(RUNNR-1)+RUNNR}=sess_wmeanBOLD{RUNNR}(vx,vy,vz,:);  %#ok<*IDISVAR,*NODEF>
        stimulus{(RUNNR-1)+RUNNR}=[];
        for voln = 1:size(stim.norm{RUNNR},2)
            stimulus{(RUNNR-1)+RUNNR} = cat(3, stimulus{(RUNNR-1)+RUNNR}, stim.norm{RUNNR}{voln}); %#ok<*AGROW>
        end
        if exist('sess_wmeanBOLD_inv') && ~isempty(sess_wmeanBOLD_inv{RUNNR}) %#ok<*EXIST>
            fmri_data{RUNNR+RUNNR}=sess_wmeanBOLD_inv{RUNNR}(vx,vy,vz,:);
            stimulus{RUNNR+RUNNR}=[];
            for voln = 1:size(stim.inv{RUNNR},2)
                stimulus{RUNNR+RUNNR} = cat(3, stimulus{RUNNR+RUNNR}, stim.inv{RUNNR}{voln});
            end
        end
    end
    
    % clear empty cell for non-existing inverse stimuli
    if isempty(fmri_data{2})
        fmri_data(2)=[]; stimulus(2)=[];
    end
    
    % copy options from results, but exclude voxel mask
    options = result.options; options.vxs=[];
    TR = 2.5;
    
    % reshape data
    for i = 1:length(fmri_data)
        data{i}=squish(fmri_data{i},4)';
    end
    
    % run analyzePRF tool for this voxel only
    mri_result = analyzePRF_modeltype(stimulus,data,TR/2,options,modeltype);
end

%% Plot fit prediction and data together ==================================
if ~LoadedResults
    res = [160 160]; resmx = max(res);
    hrf = mri_result.options.hrf; % HRF that was used in the model
    degs = mri_result.options.maxpolydeg;  % max polynomial deg
    
    % Pre-compute cache for faster execution
    [d,xx,yy] = makegaussian2d(resmx,2,2,2,2);
    
    % Prepare the stimuli for use in the model
    stimulusPP = {};
    for p=1:length(stimulus)
        stimulusPP{p} = squish(stimulus{p},2)';
        stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];
    end
    
    % Define the model function.
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
    
    % Construct projection matrices that fit and remove the polynomials.
    polymatrix = {};
    for p=1:length(degs)
        polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
    end
    
    % Which voxel should we inspect?
    vx = 1; % we only do one here, but let's stay flexible
    
    % Collect the data and the model fit. Remove slow trends in the data.
    datats = {}; modelts = {};
    for p=1:length(data)
        datats{p} =  polymatrix{p}*data{p}(vx,:)';
        modelts{p} = polymatrix{p}*modelfun(mri_result.params(1,:,vx),single(stimulusPP{p}));
    end
end

%% Save the results =======================================================
if ~LoadedResults
    save('ExampleVoxel','datats','modelts');
end

%% Visualize the results ==================================================
f=figure; hold on;
set(gcf,'Units','points','Position',[100 100 1000 300]);
plot(datats{4},'ok','MarkerSize',6,'MarkerFaceColor',[.75 .75 .75])
plot(modelts{4},'k','Linewidth',2)
%plot(cat(1,datats{:}),'r-');
%plot(cat(1,modelts{:}),'b-');
%straightline(300*(1:4)+.5,'v','g-');
xlabel('Time (stimulus position)'); ylabel('BOLD signal');
ax = axis; 
axis([.5 460+.5 ax(3:4)]);
title('Time-series data MRI - Example voxel');
legend({'Data','Model'})

saveas(f,fullfile(pngfld, 'ExampleVoxel.png'));
saveas(f,fullfile(svgfld, 'ExampleVoxel.svg'));

%close all;