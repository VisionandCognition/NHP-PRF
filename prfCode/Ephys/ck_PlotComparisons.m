function ck_PlotComparisons(subj,sess)
% plots the PRF fitting results for MUA data
% c.klink@nin.knaw.nl

%% init ===================================================================
if nargin<2
    subj = 'Lick';
    sess = '20180807_B2';
    fprintf(['Using defaults, SUBJ: ' subj ', SESS: ' sess '\n\n']);
end

%% data location ==========================================================
base_path = pwd;
if ismac % laptop & portable HDD
    data_path = '/Volumes/MRI_2/PRF_EPHYS';
    home_path = '/Users/chris';
else % 
    data_path = '/media/NETDISKS/VCNIN/PRF_EPHYS';
    home_path = '/home/chris';
end
raw_fld = fullfile(data_path,'Data_raw');
pp_fld = fullfile(data_path,'Data_preproc');
beh_fld = fullfile(data_path,'Log_pRF');
proc_fld = fullfile(data_path,'Data_proc');

%% load ===================================================================
load(fullfile(proc_fld,subj,sess,'MUA',...
    [subj '_' sess '_allarrays_PRFMUA']));
load(fullfile(proc_fld,subj,sess,'LFP',...
    [subj '_' sess '_allarrays_PRFLFP']));


%% Histograms of R2 
figure
subplot(2,3,1);hold on;box on;
rrr=[];for i=1:8; rrr=[rrr;PRF_mua(i).rcrf'];end
hist(rrr,0:0.05:10);
title('Histogram R2 size MUA')

for fb=1:5
    subplot(2,3,fb+1);hold on;box on;
    rrr=[];for i=1:8; rrr=[rrr;PRF_lfp(i,fb).rcrf'];end
    hist(rrr,0:0.05:10)
    title(['Histogram R2 size LFP: ' FreqLabel{fb}])
end

%% MUA vs LFP size
mvl_s=[];mvl_r=[];
for a=1:8
    mvl_s=[mvl_s; PRF_mua(a).rsd(1:128) PRF_lfp(a,1).rsd PRF_lfp(a,2).rsd ...
        PRF_lfp(a,3).rsd PRF_lfp(a,4).rsd PRF_lfp(a,5).rsd];
    mvl_r=[mvl_r; PRF_mua(a).rcrf(1:128)' PRF_lfp(a,1).rcrf' PRF_lfp(a,2).rcrf' ...
        PRF_lfp(a,3).rcrf' PRF_lfp(a,4).rcrf' PRF_lfp(a,5).rcrf'];    
end
figure;
for i=1:5
    subplot(2,3,i); hold on; box on;
    plot([0 6],[0 6],'k','LineWidth',1);
    inc=mvl_r(:,i+1)>0.5;
    scatter(mvl_s(inc,1),mvl_s(inc,i+1));
    set(gca,'xlim',[0 6],'ylim',[0 6],'FontSize',14);
    title(['MUA vs LFP-' FreqLabel{fb}],'FontSize',20)
end


