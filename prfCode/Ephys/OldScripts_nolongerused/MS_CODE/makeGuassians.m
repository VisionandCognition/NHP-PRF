%% Load the data & log
load('pRF_Bumba_20180105')
%Loads up the response to the 260 bars from the 192 channels (bars x chans)
% MUAm > mean response of each channel to each stimulus bar (averages over
% time window ???)

%% Create stimulus masks
%MAke stimuli matrices (260 bars)
[px,py] = meshgrid(1:1024,1:768);
px = px-512;
py = (769-py)-384;
M = zeros(size(px,1),size(px,2),260);
for n = 1:260
    xpos = LOG.barx(LOG.details(n,1));
    ypos = LOG.bary(LOG.details(n,2));
    ori = LOG.details(n,3);
    bx = LOG.bar(:,1,ori)+xpos;
    by = LOG.bar(:,2,ori)+ypos;
    M(:,:,n) = inpolygon(px,py,bx,by);
end

%% MAKE GAUSSIANS AND TIMECOURSES
lsz = 1024;
%In degrees
lszdeg = lsz/pixperdeg;
fwhm = 2*sqrt(log(2)*2);

x = ((1:1024)-512)./pixperdeg;
y = ((768:-1:1)-384)./pixperdeg;

ix = -2:0.25:8;
iy = 2:-0.25:-8;
sd = ([0.2:0.2:6]./fwhm);
na = length(ix);
nb = length(iy);
nc = length(sd);

gdets = NaN(na*nb*nc,3);
gz = 0;
for a = 1:na
    for b = 1:nb
        for c = 1:nc
            gz = gz+1;
            gdets(gz,:) = [ix(a),iy(b),sd(c)];
        end
    end
end
ngauss = size(gdets,1);

%% Multiply gaussians with stimulus masks to predict responses
%First make the Gaussians for all tehse positions
nconds = 260;
TC = zeros(nconds,ngauss);
tic
my_pp=parpool(2); %%%%% NUMBER OF CORES!!! RUN ON CLUSTER?
parfor gz = 1:10%ngauss %PAR
    
    %Read in Gaussian parameters
    a = gdets(gz,1);
    b = gdets(gz,2);
    c = gdets(gz,3);
    %MAke Gaussian
    gw = normpdf(x,a,c);
    gv = normpdf(y,b,c);
    Z = gv'*gw;
    Z = Z./sum(sum(Z)); %Could normalise by sum
    
    %Multiply Gaussian by stimulus maps to get predicted response to each
    %stimulus
    for cnd = 1:nconds
        TC(cnd,gz) = sum(sum(Z.*M(:,:,cnd)));
    end
end
toc

%% Save the predicted response timecourses for the gaussians
save('TimeCourses_Bumba_20180105','TC','gdets')