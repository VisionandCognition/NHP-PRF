function ck_PlotPRF_LFP(subj,sess)
% plots the PRF fitting results for LFP data
% split by frequency band
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
load(fullfile(proc_fld,subj,sess,'LFP',...
    [subj '_' sess '_allarrays_PRFLFP']));

%% Plot ===================================================================
r_min = 0.5;

%% visual field coverage
FreqLabel={'Theta','Alpha','Beta','Gamma1','Gamma2'};
LFP_pfit = []; 
for fb=1:5
    coll=[];
    ff1(fb)=figure;hold on;
    for a=1:8
        subplot(2,4,a); hold on;box on;
        plot([0 0],[-100 100],'k','linewidth',.5);
        plot([-100 100],[0 0],'k','linewidth',.5);
        inc = PRF_lfp(a,fb).rcrf>r_min & PRF_lfp(a,fb).rsd'~=max(PRF_lfp(a,fb).rsd);
        viscircles([PRF_lfp(a,fb).rmx(inc),PRF_lfp(a,fb).rmy(inc)],...
            PRF_lfp(a,fb).rsd(inc)./2,'Color','r','LineWidth',1);
        coll=[coll; sqrt(PRF_lfp(a,fb).rmx(inc).^2+PRF_lfp(a,fb).rmy(inc).^2) ...
            PRF_lfp(a,fb).rsd(inc) PRF_lfp(a,fb).rcrf(inc)'];
        set(gca,'xlim',[-4 12],'ylim',[-12 4],'FontSize',14);
        title(['Instance ' num2str(a) ' ' FreqLabel{fb}],'FontSize', 18); 
    end
    set(ff1(fb),'Position',[0 0 1200 600])
    
    % do a fit
    [pf1,s1] = polyfit(coll(:,1),coll(:,2),1);
    LFP_pfit = [LFP_pfit; pf1];
    save('~/Desktop/LFP_polyfit','LFP_pfit');

    % Figure of average size per eccentricity and full scatter
    b=[];
    for ecc = 1:8
        sel = coll(:,1)>ecc-1 & coll(:,1)<ecc;
        b=[b; ecc-.5 mean(coll(sel,2)) std(coll(sel,2))./sqrt(sum(sel))];
    end
    figure; subplot(1,2,1);
    errorbar(b(:,1),b(:,2),b(:,3),'o');
    title(['Mean pRF size / Ecc. ,' FreqLabel{fb}]);

    subplot(1,2,2);
    scatter(coll(:,1),coll(:,2));
    title(['pRF size vs Ecc, ' FreqLabel{fb}]);
end

%% ECC vs Size
for fb=1:5
    %
    ff(fb)=figure;hold on;
    fff(fb)=figure;hold on;
    for a=1:8
        figure(ff(fb));
        subplot(2,4,a);hold on;box on; 
        inc = PRF_lfp(a,fb).rcrf>r_min & PRF_lfp(a,fb).rsd'~=max(coll(:,2));
        ecc = sqrt(PRF_lfp(a,fb).rmx(inc).^2+PRF_lfp(a,fb).rmy(inc).^2);
        rfsz = PRF_lfp(a,fb).rsd(inc);
        plot(ecc,rfsz,'o','MarkerSize',8)
        set(gca,'FontSize',14);
        title(['ECC vs Size, Inst. ' num2str(a) ', ' FreqLabel{fb}],'FontSize', 18);
        set(ff(fb),'Position',[0 0 1200 600])

        figure(fff(fb))
        plot(ecc,rfsz,'o','MarkerSize',8);
        box on;
        set(gca,'FontSize',14,'xlim',[0 12],'ylim',[0 8])
        title(['Ecc vs Size ALL arrays ,' FreqLabel{fb}])
        
    end
end

%% Example Arrays
exarr=6;
exchan=1;
figure; hold on;

for fb=1:5
    subplot(5,1,fb);hold on;box on;
    bar(PRF_lfp(exarr,fb).LFPm(:,exchan)-min(PRF_lfp(exarr,fb).LFPm(:,exchan)),...
        'BarWidth',1,'FaceColor','k')
    for sw=0:30:240
        plot([sw sw],[-0.5 3],'Color',[0 .5 .75],'LineWidth',2);
    end
    %set(gca,'xlim',[0 240],'ylim',[0 1.1],'FontSize',14);
    if fb==1
        title(['LFP response / bar [Example channel] ' FreqLabel{fb}],'FontSize',20)
    elseif fb==5
        xlabel('Bar position','FontSize',14);
        title(FreqLabel{fb},'FontSize',20);
    end
end

