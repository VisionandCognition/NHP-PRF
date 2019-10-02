x=55; y=15; z=45;

idx=find(result.R2==max(max(max(result.R2))));

tseries=[];
for t=1:length(s_run(3).vol)
    tseries=[tseries s_run(3).vol{t}(idx)];
end

% tseries is BOLD timeseries



%% concatenate -----
stimulus={};fmri_data={};
doExtraRegression=true;
fprintf('Concatenating stimuli and volumes...\n');
for r=1:length(s_run)
    stimulus{r}=[]; fmri_data{r}=[];
    for voln = 1:size(s_run(r).vol,2)
        stimulus{r} = cat(3, stimulus{r}, s_run(r).stim{voln}); %#ok<*AGROW>
        fmri_data{r} = cat(4, fmri_data{r}, s_run(r).vol{voln});
    end
    if doExtraRegression
        extraregr{r} = cat(2,s_run(r).motion.estimates,s_run(r).rew);
        % 12 motion correction parameters + reward events
    end
end


% get the prediction
% Define some variables
res = [160 160];                    % row x column resolution of the stimuli
resmx = 160;                        % maximum resolution (along any dimension)
hrf = result.options.hrf;          % HRF that was used in the model
degs = result.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);
stimulusPP = {};
for p=1:length(stimulus)
  stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

for d=1:length(fmri_data)
    data{d}=[];
    for d2=1:size(fmri_data{d},4)
        data{d}=[data{d} reshape(fmri_data{d}(:,:,:,d2),numel(fmri_data{d}(:,:,:,d2)),1)];
    end
end


%data=fmri_data;

polymatrix = {};
for p=1:length(degs)
  polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},2),0:degs(p)));
end

vx=idx;


% For each run, collect the data and the model fit.  We project out polynomials
% from both the data and the model fit.  This deals with the problem of
% slow trends in the data.
datats = {};
modelts = {};
for p=1:length(data)
  datats{p} =  data{p}(vx,:)';
  modelts{p} = modelfun(result.params(1,:,vx),stimulusPP{p});
end

% Visualize the results
figure; hold on;
set(gcf,'Units','points','Position',[100 100 1000 100]);
plot(cat(1,datats{:}),'r-');
plot(cat(1,modelts{:}),'b-');
straightline(300*(1:4)+.5,'v','g-');
xlabel('Time (s)');
ylabel('BOLD signal');
ax = axis;
axis([.5 1200+.5 ax(3:4)]);
title('Time-series data');






%%
TR=2.5;
options.wantglmdenoise = extraregr;
result = analyzePRF(stimulus,data2,TR/2,options);

cfactor = 10/100;




