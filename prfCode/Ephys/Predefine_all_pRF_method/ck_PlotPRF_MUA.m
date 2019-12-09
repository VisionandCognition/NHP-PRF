function ck_PlotPRF_MUA(subj,sess)
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
    data_path = '/media/NETDISKS/VS02/VandC/PRF_EPHYS';
    home_path = '/home/chris';
end
raw_fld = fullfile(data_path,'Data_raw');
pp_fld = fullfile(data_path,'Data_preproc');
beh_fld = fullfile(data_path,'Log_pRF');
proc_fld = fullfile(data_path,'Data_proc');

%% load ===================================================================
load(fullfile(proc_fld,subj,sess,'MUA',...
    [subj '_' sess '_allarrays_PRFMUA']));

%% Plot ===================================================================
r_min = 0.50;

%% visual field coverage
f1=figure;hold on;
coll=[];
for a=[1:8]
    subplot(2,4,a); hold on; box on;
    plot([0 0],[-100 100],'k','linewidth',.5);
    plot([-100 100],[0 0],'k','linewidth',.5);
    inc = PRF_mua(a).rcrf>r_min & PRF_mua(a).rsd'~=max(PRF_mua(a).rsd);
    viscircles([PRF_mua(a).rmx(inc),PRF_mua(a).rmy(inc)],...
        PRF_mua(a).rsd(inc)./2,'Color','r','LineWidth',1);
    coll=[coll; sqrt(PRF_mua(a).rmx(inc).^2+PRF_mua(a).rmy(inc).^2) ...
        PRF_mua(a).rsd(inc) PRF_mua(a).rcrf(inc)'];
    set(gca,'xlim',[-4 12],'ylim',[-12 4],'FontSize',14);
    title(['Instance ' num2str(a)],'FontSize', 18); 
end
set(f1,'Position',[0 0 1200 600])

%% ECC vs Size
f2=figure;hold on;
f3=figure;hold on;
f4=figure;hold on;
coll=[];
for a=[1:8]
    figure(f2)
    subplot(2,4,a); hold on; box on;
    inc = PRF_mua(a).rcrf>r_min & PRF_mua(a).rsd'~=max(PRF_mua(a).rsd);
    ecc = sqrt(PRF_mua(a).rmx(inc).^2+PRF_mua(a).rmy(inc).^2);
    rfsz = PRF_mua(a).rsd(inc);
    plot(ecc,rfsz,'o','MarkerSize',8)
    set(gca,'FontSize',14);
    title(['ECC vs Size, Inst. ' num2str(a)],'FontSize', 18);
    set(f2,'Position',[0 0 1200 600])
    
    figure(f3)
    plot(ecc,rfsz,'o','MarkerSize',8)
    
    figure(f4);
    subplot(2,4,a); hold on; box on;
    hist(rfsz./ecc,0:.1:6);

    coll=[coll;rfsz ecc];
end
[pf1,s1] = polyfit(coll(:,2),coll(:,1),1);
MUA_pfit = pf1;
save('~/Desktop/MUA_polyfit','MUA_pfit');
    
%% histogram of slope
figure;hist(coll(:,1)./coll(:,2),0:.1:2)
title('Histogram of SIZE/ECC')

%% visual field coverage all channels
figure;hold on;
for a=[1:8]
    %subplot(2,4,a); 
    hold on;box on;
    plot([0 0],[-100 100],'k','linewidth',.5);
    plot([-100 100],[0 0],'k','linewidth',.5);
    inc = PRF_mua(a).rcrf>r_min & PRF_mua(a).rsd'~=max(coll(:,2));
    plot(PRF_mua(a).rmx(inc),PRF_mua(a).rmy(inc),'*');
    set(gca,'xlim',[-3 9],'ylim',[-9 3]);
    title('Total visual field coverage','FontSize',18)
    set(gca,'FontSize',14);
end

% b=[];
% for ecc = 1:8
%   sel = coll(:,1)>ecc-1 & coll(:,1)<ecc;
%   b=[b; ecc-.5 mean(coll(sel,2)) std(coll(sel,2))./sqrt(sum(sel))];
% end
% figure; subplot(1,2,1);
% errorbar(b(:,1),b(:,2),b(:,3),'o')
% 
% subplot(1,2,2);
% scatter(coll(:,1),coll(:,2))

%% example of responses for bar sweeps
exarr=1;
exchan=1;
figure; hold on; box on;
bar(PRF_mua(exarr).MUAm(:,exchan)-min(PRF_mua(exarr).MUAm(:,exchan)),...
    'BarWidth',1,'FaceColor','k')
for sw=0:30:240
    plot([sw sw],[-0.5 3],'Color',[0 .5 .75],'LineWidth',2);
end
title('MUA response per bar position [Example channel]','FontSize',20)
set(gca,'xlim',[0 240],'ylim',[0 1.1],'FontSize',14);
xlabel('Bar position','FontSize',14);
ylabel('MUA (mean of 5 runs)','FontSize',14);

%% histogram of R2 values
figure
hold on;box on;
rrr=[];for i=1:8; rrr=[rrr;PRF_mua(i).rcrf'];end
hist(rrr,0:0.05:1)
