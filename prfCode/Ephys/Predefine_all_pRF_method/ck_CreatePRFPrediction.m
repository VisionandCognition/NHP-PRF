function ck_CreatePRFPrediction(subj,sess)
% Create PRF predictions
% multiply stimuli with a bunch of gaussians
% find best fitting gaussian using ck_FitPRF.m
% c.klink@nin.knaw.nl

% ------ SHOULDN'T MATTER BUT WE CAN BASE THIS ON DIFFERENT LOGS -----
if nargin<2
    subj = 'Lick';
    sess = '20180807_B2';
    fprintf(['Using defaults, SUBJ: ' subj ', SESS: ' sess '\n\n']);
end
% --------------------------------------------------------------------

fprintf('=============================================================\n');
fprintf(' Loading and processing stimuli for %s, %s\n', subj, sess);
fprintf('=============================================================\n');

%% data location
base_path = pwd;
if ismac % laptop & portable HDD
    data_path = '/Volumes/MRI_2/PRF_EPHYS';
elseif ispc
    data_path = '\\vcnin\PRF_EPHYS';
else % 
    data_path = '/media/NETDISKS/VS02/VandC/PRF_EPHYS';
end
raw_fld = fullfile(data_path,'Data_raw');
pp_fld = fullfile(data_path,'Data_preproc');
beh_fld = fullfile(data_path,'Log_pRF');
save_fld = fullfile(data_path,'Data_proc');

%% Load stimulus
cd(fullfile(save_fld,subj,sess));
warning off
load([subj '_' sess '_STIM']);
warning on
cd(base_path);
SM = cell2mat(STIM.img);
SM = reshape(SM,[size(STIM.img{1}) length(STIM.img)]);
cd(base_path);

%% Create Gaussians
fprintf('Creating Gaussians as predictions for RFs...\n');
lsz = size(STIM.img{1},1);

%In degrees
lszdeg = lsz/B(1).Par.PixPerDeg;
fwhm = 2*sqrt(log(2)*2);

x = ((1:lsz)-lsz/2)./B(1).Par.PixPerDeg;
y = ((lsz:-1:1)-lsz/2)./B(1).Par.PixPerDeg;

radius = 15;
% Original >> huge range >> takes long
%ix = -radius:0.2:radius;
%iy = radius:-0.2:-radius;

% Reduced range based on earlier results and anatomy
%ix = -3:0.2:radius;
%iy = 3:-0.2:-radius;
%sd = 0.1:0.1:8;
ix = -1.5:0.1:10;
iy = 1:-0.1:-12;
sd = 0.1:0.1:6;

fprintf(['VF Coverage for prf predictions is X: ' num2str(min(ix)) ' to ' ...
     num2str(max(ix)) ', Y: ' num2str(min(iy)) ' to ' num2str(max(iy)) '\n']);
fprintf(['RF size Coverage for prf predictions is: ' num2str(min(sd)) ' to ' ...
     num2str(max(sd)) '\n']);

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
fprintf('Number of Gaussians created: %d\n', ngauss);

%% Predict gaussian responses 
fprintf('Multiplying Gaussian with stim-masks (will take a while)\n')

% NB! This is computationally intensive
ncores = 12; % do this thing in parallel on how many cores?
fprintf('Will do this in parallel on %d cores\n', ncores')

% First make the Gaussians for all positions
nconds = size(SM,3);
PRF_pred = zeros(nconds,ngauss);
delete(gcp('nocreate'))
my_pp=parpool(ncores); %%%%% NUMBER OF CORES!!! RUN ON CLUSTER FOR SPEED

%parfor gz = 1:10 % test with only a few
parfor gz = 1:ngauss %PAR
    tic
    %Read in Gaussian parameters
    a = gdets(gz,1);
    b = gdets(gz,2);
    c = gdets(gz,3);
    %Make Gaussian
    gw = normpdf(x,a,c);
    gv = normpdf(y,b,c);
    Z = gv'*gw;
    Z = Z./sum(sum(Z)); %Could normalise by sum
    
    Z2=reshape(repmat(Z,[1,nconds]),[size(Z) nconds]);
    
    toc
    %Multiply Gaussian by stimulus maps to get predicted response to each
    %stimulus
    for cnd = 1:nconds
        PRF_pred(cnd,gz) = sum(sum(Z.*double(SM(:,:,cnd))));
    end
    toc
    fprintf(['Created prediction number ' num2str(gz) '\n']);
end
fprintf('\nDone!\n');
delete(gcp('nocreate'));

%% save predicted responses
fprintf('Saving the predicted PRF responses...\n');
save(fullfile(data_path,'PRF_Prediction'),'PRF_pred','SM','gdets','-v7.3');