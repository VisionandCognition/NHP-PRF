load('pRF_Bumba_20180105')
load('TimeCourses_Bumba_20180105')
nchans = 192;
for s = 1:12
    slice(s).ind = ones(nchans,1);
    slice(s).bestbeta = NaN(2,nchans);
    slice(s).currentbest = inf(1,nchans);
end

%Divide the timecourse up into slices
ngauss = size(TC,2);
nslices = 12;
d = floor(linspace(1,ngauss,nslices+1));
tic
parfor s = 1:nslices
    st = d(s);
    ed = d(s+1)-1;
    for gz = st:1:ed
        %Make design matrix for this predicted timecourse (no constant
        %term!)
        %Fit the GLM using the matrix form
        %b = pinv(dm)*pixmat; %pinv = inv(dm'*dm)*dm'
        dm = [TC(:,gz),ones(260,1)];
        [b,~,mss] = lscov(dm,MUAm); %get best fitting betas
        mss(b(1,:)<0) = Inf; %Set negative betas to Inf
        %Compare sse to current best, note down impriovements in teh minimum
        slice(s).ind(mss<slice(s).currentbest) = gz;
        slice(s).bestbeta(:,mss<slice(s).currentbest) = b(:,mss<slice(s).currentbest);
        slice(s).currentbest(mss<slice(s).currentbest) = mss(mss<slice(s).currentbest); %Update currentbest
    end
end
toc

%Combine slices
best = zeros(nslices,nchans);
allbeta = zeros(2,nchans,nslices);
allind = zeros(nslices,nchans);
for s = 1:nslices
    best(s,:) = slice(s).currentbest;
    allbeta(:,:,s) = slice(s).bestbeta;
    allind(s,:) = slice(s).ind;
end
[bestval,bestix] = min(best);
clear bestbeta
clear ind
for i = 1:nchans
    ind(i) = allind(bestix(i),i);
    bestbeta(:,i) = allbeta(:,i,bestix(i));
end

%READ in values
rmx = gdets(ind,1);
rmy = gdets(ind,2);
rsd = gdets(ind,3);

%Convert sd back into FWHM
rsd_fwhm = rsd.*fwhm;

 %Calculate corrcoef
 %turn into matrix computation
 rcrf = NaN(1,nchans);
 for N = 1:nchans
     %Get this channels
     pixel = MUAm(:,N);
     dm = [TC(:,ind(N)),ones(260,1)]; %No offset allowed!
     %Fit the GLM using the matrix form
     b = pinv(dm)*pixel; %pinv = inv(dm'*dm)*dm'
     %Gives a best-fitting beta for each pixel.
     %Make the predicted response using the GLM equation y = x*beta
     y = dm*b;
     c = corrcoef(y,MUAm(:,N));
     rcrf(N) = c(1,2);
 end
 
 save('Bumba_pRF2','rmx','rmy','rsd','rcrf')
 
 if 1
    %Example fit
    [i,j] = max(rcrf);
    %Get response from this pixel
    bestpix = MUAm(:,j);
    bb = bestbeta(:,j);
    %Get the predicted time-course for the best-fitting Gaussian
    %REMEMBER ind is the index of the best-fitting timecourse for each pixel
    predTC = TC(:,ind(j)).*bb(1)+bb(2);
    figure,bar(1:260,bestpix),hold on,scatter(1:260,predTC,'r','filled')
 end

 figure
 for j = 1:nchans
  bestpix = MUAm(:,j);
    bb = bestbeta(:,j);
    %Get the predicted time-course for the best-fitting Gaussian
    %REMEMBER ind is the index of the best-fitting timecourse for each pixel
    predTC = TC(:,ind(j)).*bb(1)+bb(2);
    bar(1:260,bestpix),hold on,scatter(1:260,predTC,'r','filled')
    ylim([bb(2)-bb(1).*1.5 bb(2)+bb(1).*1.5])
    drawnow
    pause
    hold off
 end
