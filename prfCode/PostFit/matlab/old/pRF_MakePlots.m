MONKEY='danny';

Thr = 5; % the R2 threshold per session
pVox = 0.2; % inclusion of voxels over sessions

%% Add nifti reading toolbox
tool_basepath = '~/Dropbox/MATLAB_NONGIT/TOOLBOX';
addpath(genpath(fullfile(tool_basepath, 'NIfTI')));

%% Define paths
respath = ['~/Documents/MRI_Analysis/pRF-NHP/Results/' MONKEY];
% respath = ['~/Documents/MRI_Analysis/pRF-NHP/Results/LISA/' MONKEY];

if strcmp(MONKEY,'danny')
    roipath = ['~/Documents/MRI_Analysis/pRF-NHP/Results/Reference/' ...
        'danny/output_files/ROI_adj/1mm/nii'];
elseif strcmp(MONKEY,'eddy')
    roipath = ['~/Documents/MRI_Analysis/pRF-NHP/Results/Reference/' ...
        'eddy/output_files/ROI/1mm/nii'];
else
    fprintf(['Unknown monkey name ' MONKEY '\n']);
end

%% load results
%load(fullfile(respath,'AllFitResults_us_motreg'));
ECC = load_nii(fullfile(respath,'AveragedResults',['Thr_' num2str(Thr)], ...
    'nii',['MeanEccentricity_th' num2str(Thr) '_pv' num2str(pVox) '.nii']));
ANG = load_nii(fullfile(respath,'AveragedResults',['Thr_' num2str(Thr)], ...
    'nii',['MeanAngle_th' num2str(Thr) '_pv' num2str(pVox) '.nii']));
RFS = load_nii(fullfile(respath,'AveragedResults',['Thr_' num2str(Thr)], ...
    'nii',['MeanRFS_th' num2str(Thr) '_pv' num2str(pVox) '.nii']));

ecc = reshape(ECC.img,[1 numel(ECC.img)]);
ang = reshape(ANG.img,[1 numel(ANG.img)]);
rfs = reshape(RFS.img,[1 numel(RFS.img)]);

ecc(ecc==-99)=nan;
ang(ang==-99)=nan;
rfs(rfs==-99)=nan;
[rfx,rfy]=pol2cart(ang*(pi/180),ecc);

%% overall plots
% scatter(ecc,rfs) ===
% bin it
ed=ecc; rfd=rfs;
rfd(isinf(rfd))=nan;
m=[]; sd=[]; sem=[];
bw=0.5; bins = .25:bw:9.25;
for i= bins
    m=[m nanmean(rfd(ed>i-bw/2 & ed<=i+bw/2))];
    sd=[sd nanstd(rfd(ed>i-bw/2 & ed<=i+bw/2))];
    sem=[sem nanstd(rfd(ed>i-bw/2 & ed<=i+bw/2))/...
        sqrt(sum(ed>i-bw/2 & ed<=i+bw/2))];
end

f=figure;
%errorbar(bins,m,sem)
bar(bins,m,'BarWidth',.9,'FaceColor',[0 0.45 0.74])
%plot(bins,sem)
xlabel('Eccentricity (dva)','FontSize',16)
ylabel('Mean pRF size','FontSize',16)
title('All Voxels: Ecc v. Size','FontSize',18)
set(gca,'xlim',[0 9.5],'FontSize',14);

%[pf1,s1] = polyfit(bins,m,1);
%hold on;
%plot([0.5 10],pf1(2)+pf1(1)*[0.5 10])

set(f,'Position',[0 0 700 900])

%% load ROIs =====
V1 = load_nii(fullfile(roipath,'V1_roi.nii'));
V2 = load_nii(fullfile(roipath,'V2_roi.nii'));
V3A = load_nii(fullfile(roipath,'V3A_roi.nii'));
V3d = load_nii(fullfile(roipath,'V3d_roi.nii'));
V3v = load_nii(fullfile(roipath,'V3v_roi.nii'));
V4 = load_nii(fullfile(roipath,'V4_roi.nii'));
V4t = load_nii(fullfile(roipath,'V4t_roi.nii'));
V4v = load_nii(fullfile(roipath,'V4v_roi.nii'));
MST = load_nii(fullfile(roipath,'MST_roi.nii'));
MT = load_nii(fullfile(roipath,'MT_roi.nii'));
TEO = load_nii(fullfile(roipath,'TEO_roi.nii'));
TPO = load_nii(fullfile(roipath,'TPO_roi.nii'));
VIP = load_nii(fullfile(roipath,'VIP_roi.nii'));
area7A = load_nii(fullfile(roipath,'7a_(Opt-PG)_roi.nii'));
area5 = load_nii(fullfile(roipath,'5_(PE)_roi.nii'));
area5a = load_nii(fullfile(roipath,'5_(PEa)_roi.nii'));
LIPd = load_nii(fullfile(roipath,'LIPd_roi.nii'));
LIPv = load_nii(fullfile(roipath,'LIPv_roi.nii'));

ROI(1).name = 'V1';
ROI(1).vox = reshape(V1.img,[1 numel(V1.img)]);
ROI(2).name = 'V2';
ROI(2).vox = reshape(V2.img,[1 numel(V2.img)]);
ROI(3).name = 'V3A';
ROI(3).vox = reshape(V3A.img,[1 numel(V3A.img)]);
ROI(4).name = 'V3d';
ROI(4).vox = reshape(V3d.img,[1 numel(V3d.img)]);
ROI(5).name = 'V3v';
ROI(5).vox = reshape(V3v.img,[1 numel(V3v.img)]);
ROI(6).name = 'V4';
ROI(6).vox = reshape(V4.img,[1 numel(V4.img)]);
ROI(7).name = 'V4t';
ROI(7).vox = reshape(V4t.img,[1 numel(V4t.img)]);
ROI(8).name = 'V4v';
ROI(8).vox = reshape(V4v.img,[1 numel(V4v.img)]);
ROI(9).name = 'MST';
ROI(9).vox = reshape(MST.img,[1 numel(MST.img)]);
ROI(10).name = 'MT';
ROI(10).vox = reshape(MT.img,[1 numel(MT.img)]);
ROI(11).name = 'TEO';
ROI(11).vox = reshape(TEO.img,[1 numel(TEO.img)]);
ROI(12).name = 'TPO';
ROI(12).vox = reshape(TPO.img,[1 numel(TPO.img)]);
ROI(13).name = 'VIP';
ROI(13).vox = reshape(VIP.img,[1 numel(VIP.img)]);
ROI(14).name = 'area7A';
ROI(14).vox = reshape(area7A.img,[1 numel(area7A.img)]);
ROI(15).name = 'area5';
ROI(15).vox = reshape(area5.img,[1 numel(area5.img)]);
ROI(16).name = 'area5a';
ROI(16).vox = reshape(area5a.img,[1 numel(area5a.img)]);
ROI(17).name = 'LIPd';
ROI(17).vox = reshape(LIPd.img,[1 numel(LIPd.img)]);
ROI(18).name = 'LIPv';
ROI(18).vox = reshape(LIPv.img,[1 numel(LIPv.img)]);

%% Plot ECC-SZ results per ROI =====
% f1=figure; 
% subplot(1,2,1);hold on;
% subplot(1,2,2);hold on;
% 
% bw=1; bins = 0.5:bw:8.5;
% for i=1:length(ROI)
%     LEG{i}=ROI(i).name;
%     % bin it
%     ed=ecc(logical(ROI(i).vox)); 
%     rfd=rfs(logical(ROI(i).vox));
%     rfd(isinf(rfd))=nan;
%     m=[];sd=[];sem=[];
%     for ii= bins
%         %m=[m nanmean(rfd(ed>ii-bw/2 & ed<=ii+bw/2))];
%         m=[m median(rfd(ed>ii-bw/2 & ed<=ii+bw/2))];
%         sd=[sd nanstd(rfd(ed>ii-bw/2 & ed<=ii+bw/2))];
%         sem=[sem nanstd(rfd(ed>ii-bw/2 & ed<=ii+bw/2))/...
%             sqrt(sum(ed>ii-bw/2 & ed<=ii+bw/2))];
%     end
%     subplot(1,2,1);
%     %errorbar(bins,m,sem);
%     plot(bins,m);
%     [pf1,s1] = polyfit(bins(1:5),m(1:5),1);
%     subplot(1,2,2);
%     plot([0 9],pf1(2)+pf1(1)*[0 9])
% end
% legend(LEG);

%% Plot ECC-SZ fit results per ROI =====
f1=figure; 
subplot(1,2,1);hold on;
subplot(1,2,2);hold on;

bw=1; bins = 0.5:bw:8.5; l=1;
for i=1:length(ROI)
    LEG{l}=ROI(i).name;l=l+1;
    % bin it
    ed=ecc(logical(ROI(i).vox)); 
    rfd=rfs(logical(ROI(i).vox));
    
    if strcmp(ROI(i).name,'V1')
        [pf1,s1] = polyfit(bins(1:5),m(1:5),1);
        BOLD_pfit = pf1;
        save(['~/Desktop/BOLD_polyfit_' MONKEY],'BOLD_pfit')
    end
    
    percentile_cutoff = 95;
    perc_sel = rfd<=prctile(rfd,percentile_cutoff);
    
    ed2=ed(perc_sel);
    rfd2=rfd(perc_sel);
    
    m=[];sd=[];sem=[];
     
    for ii= bins
        %m=[m nanmean(rfd2(ed2>ii-bw/2 & ed2<=ii+bw/2))];
        m=[m median(rfd2(ed2>ii-bw/2 & ed2<=ii+bw/2))];
        sd=[sd nanstd(rfd2(ed2>ii-bw/2 & ed2<=ii+bw/2))];
        sem=[sem nanstd(rfd2(ed2>ii-bw/2 & ed2<=ii+bw/2))/...
            sqrt(sum(ed2>ii-bw/2 & ed2<=ii+bw/2))];
    end
    subplot(1,2,1);
    %errorbar(bins,m,sem);
    plot(bins,m,'o','MarkerSize',8,'LineWidth',2);
    set(gca,'FontSize',14)
    xlabel('Eccentricity (dva)','FontSize',16)
    title('Ecc. vs Size per ROI','FontSize',16)
    [pf1,s1] = polyfit(bins(1:5),m(1:5),1);
    %figure(f2);
    subplot(1,2,2);
    plot([0 9],pf1(2)+pf1(1)*[0 9],'LineWidth',2)
    set(gca,'FontSize',14)
    xlabel('Eccentricity (dva)','FontSize',16)
    title('Linear Fit','FontSize',16)
    %LEG{l}=[ ROI(i).name ' fit'];l=l+1;
end
subplot(1,2,1);legend(LEG,'Location','NorthWest','FontSize',14);
subplot(1,2,2);legend(LEG,'Location','NorthWest','FontSize',14);
set(f1,'Position',[0 0 1600 700])

%% Plot VISUAL FIELD coverage ============
close all; warning off;
for i=1:length(ROI)
    f3=figure; scrsize=18;
    
    % bin it
    rfxs=rfx(logical(ROI(i).vox)); 
    rfys=rfy(logical(ROI(i).vox)); 
    rfss=rfs(logical(ROI(i).vox));
    % only include data <=nth percentile in size
    percentile_cutoff = 95;
    perc_sel = rfss<=prctile(rfss,percentile_cutoff);
    
    subplot(5,2,1:2:7);hold on;
    plot(rfxs(perc_sel)',rfys(perc_sel)','o','MarkerSize',5);
    viscircles([0 0],8,'Color','k','LineWidth',3);
    t=title(ROI(i).name); set(t,'FontSize',18)
    axrange = [ get(gca,'XLim') get(gca,'YLim')];
    %set(gca,'Xlim',[-ceil(max(abs(axrange))) ceil(max(abs(axrange)))], ... 
    %    'Ylim',[-ceil(max(abs(axrange))) ceil(max(abs(axrange)))],'Box','on');
    set(gca,'Xlim',[-scrsize scrsize],'Ylim',...
        [-scrsize scrsize],'Box','on','FontSize',14);
    xlabel('H-ECC','FontSize',16); ylabel('V-ECC','FontSize',16);
    
    subplot(5,2,2:2:8);hold on;
    viscircles([rfxs(perc_sel)',rfys(perc_sel)'],...
        rfss(perc_sel)./2,'LineWidth',1);
    viscircles([0 0],8,'Color','k','LineWidth',3);
    t=title(ROI(i).name);
    set(t,'FontSize',18)
    axrange = [ get(gca,'XLim') get(gca,'YLim')];
    %set(gca,'Xlim',[-ceil(max(abs(axrange))) ceil(max(abs(axrange)))], ... 
    %    'Ylim',[-ceil(max(abs(axrange))) ceil(max(abs(axrange)))],'Box','on');
    set(gca,'Xlim',[-scrsize scrsize],'Ylim',...
        [-scrsize scrsize],'Box','on','FontSize',14);
    %xlabel('H-ECC','FontSize',16); 
    ylabel('V-ECC','FontSize',16);

    subplot(5,2,10);hold on;
    rsz = rfss(perc_sel);
    rx = rfxs(perc_sel);
    bw=2;
    xs=[];
    for x_pos=-20:bw:20
        xs=[xs; x_pos nanmean(rsz(rx>x_pos-bw/2 & rx<x_pos+bw/2))];
    end
    b=bar(xs(:,1),xs(:,2));
    set(b,'BarWidth',1,'FaceColor',[0 0.45 0.74])
    %set(gca,'Xlim',[-ceil(max(abs(axrange))) ceil(max(abs(axrange)))],'Box','on')
    set(gca,'Xlim',[-scrsize scrsize],'Box','on','FontSize',14)
    xlabel('H-ECC','FontSize',16); ylabel('pRF-size','FontSize',16);

    % save
    set(f3,'Position',[0 0 1600 800])
    savefld = fullfile(respath,'AveragedResults',['Thr_' num2str(Thr)], ...
        'VisualFieldCoverage');
    mkdir(savefld);
    subfld=fullfile(savefld,['Th' num2str(Thr) '_pv' num2str(pVox)]);
    mkdir(subfld);
    saveas(f3,fullfile(subfld,[ROI(i).name '.png']));
end
close all; warning on;

