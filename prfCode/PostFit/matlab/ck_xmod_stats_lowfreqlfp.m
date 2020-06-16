%% stats for linear models comparison

for gainopt = -1%[-1 1] % do separate for negative and positive gain
    fprintf('===================================================\n');
    fprintf(['GAIN MODE (low freq LFP): ' num2str(gainopt) '\n']);
    fprintf('===================================================\n');
    
    %% DATA ===================================================================
    clear stats
    ecc=[]; sz=[]; % eccentricity and size
    signal = []; % categorical 1=MRI, 2=MUA, 3-7=LFP bands
    area=[]; % categorical 1=V1,4=V4
    model=[]; % categorical 1=lin, 2=lin_ng, 3=css, 4=dog
    gain=[];
    
    % v1
    ecc = [ecc; ...
        mri1(m).ECC; ...
        mua1(m).ECC; ...
        ];
    sz = [ sz; ...
        mri1(m).S; ...
        mua1(m).S; ...
        ];
    signal = [signal;...
        1*ones(size(mri1(m).S)); ...
        2*ones(size(mua1(m).S)); ...
        ];
    area = [area; ...
        1*ones(size(mri1(m).S)); ...
        1*ones(size(mua1(m).S)); ...
        ];
    gain = [gain; ...
        mri1(m).G; ...
        mua1(m).G; ...
        ];
    
    for i=1:5
        if i<=3
            gainselect = ...
                lfp1(i,m).G./abs(lfp1(i,m).G) == gainopt;
            ecc = [ecc; lfp1(i,m).ECC(gainselect)];
            sz = [ sz; lfp1(i,m).S(gainselect)];
            gain = [ gain; lfp1(i,m).G(gainselect)];
            signal = [signal; (2+i)*ones(size(lfp1(i,m).ECC(gainselect)))];
            area = [area; 1*ones(size(lfp1(i,m).ECC(gainselect)))];
        else
            ecc = [ecc; lfp1(i,m).ECC];
            sz = [ sz; lfp1(i,m).S];
            gain = [ gain; lfp1(i,m).G];
            signal = [signal; (2+i)*ones(size(lfp1(i,m).S))];
            area = [area; 1*ones(size(lfp1(i,m).S))];
        end
    end
    
    % v4
    ecc = [ecc; ...
        mri4(m).ECC; ...
        mua4(m).ECC; ...
        ];
    sz = [ sz; ...
        mri4(m).S; ...
        mua4(m).S; ...
        ];
    signal = [signal;...
        1*ones(size(mri4(m).S)); ...
        2*ones(size(mua4(m).S)); ...
        ];
    area = [area; ...
        4*ones(size(mri4(m).S)); ...
        4*ones(size(mua4(m).S)); ...
        ];
    gain = [gain; ...
        mri4(m).G; ...
        mua4(m).G; ...
        ];
    
    for i=1:5
        if i<=3
            gainselect = ...
                lfp4(i,m).G./abs(lfp4(i,m).G) == gainopt;
            ecc = [ecc; lfp4(i,m).ECC(gainselect)];
            sz = [ sz; lfp4(i,m).S(gainselect)];
            gain = [ gain; lfp4(i,m).G(gainselect)];
            signal = [signal; (2+i)*ones(size(lfp4(i,m).ECC(gainselect)))];
            area = [area; 4*ones(size(lfp4(i,m).ECC(gainselect)))];
        else
            ecc = [ecc; lfp4(i,m).ECC];
            sz = [ sz; lfp4(i,m).S];
            gain = [ gain; lfp4(i,m).G];
            signal = [signal; (2+i)*ones(size(lfp4(i,m).S))];
            area = [area; 4*ones(size(lfp4(i,m).S))];
        end
    end
    
    %
    if gainopt > 0
        sigtype = {1,2,3,4,5,6,7;'mri' ,'mua','alpha-pos','beta-pos','theta-pos','hgam','lgam'};
    else
        sigtype = {1,2,3,4,5,6,7;'mri' ,'mua','alpha-neg','beta-neg','theta-neg','hgam','lgam'};
    end
    signal_num=signal;
    for i=1:length(signal_num)
        sn{i,1}=sigtype{2,signal_num(i)};
    end
    signal=categorical(signal);
    area=categorical(area);
    areaname = {' V1', 'V4'};
    DT=table(signal,area, ecc, sz, gain);
    
    %% PLOT SCATTER AND FITS ==================================================
    f_xmod = figure;
    set(f_xmod,'Position',[100 100 2000 800]);
    CLR=[...
        0 0 0;...
        0.85,0.33,0.10;...
        0.93,0.69,0.13;...
        0.47,0.67,0.19;...
        0.49,0.18,0.56;...
        0.30,0.75,0.93;...
        0.00,0.45,0.74;...
        ];
    
    % do the modelfits and plot scatter and fitlines
    % v1 -----
    DT1=DT(DT.area==categorical(1) & DT.ecc<=MaxECC,:);
    mdl_v1 = fitlm(DT1,'sz ~ 1 + signal*ecc');
    slidx = 1+ length(mdl_v1.Coefficients.Estimate)/2;
    
    v1idx=[];
    for j=1:length(mdl_v1.CoefficientNames)
        if mdl_v1.CoefficientNames{j}(1:3) == 'sig' & ...
                mdl_v1.CoefficientNames{j}(end-2:end) ~= 'ecc'
            v1idx = [v1idx str2double(mdl_v1.CoefficientNames{j}(8:end))];
        end
    end
    
    msz = 40;
    subplot(2,2,1);
    ll={}; US=[sigtype{1,:}]; hold on; ii=1;
    for i=1:length(US)
        sel=DT1.signal==categorical(i);
        if sum(sel)>0
            scatter(DT1.ecc(sel),DT1.sz(sel),msz,...
                'Marker','o',...
                'MarkerEdgeColor',CLR(i,:),'MarkerFaceColor',CLR(i,:),...
                'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.2);
            ll{ii}=sigtype{2,i};
            ii=ii+1;
        end
    end
    legend(ll,'Location','NorthEastOutside');
    title('V1');xlabel('Ecc');ylabel('Sz');
    set(gca,'xlim',[0 12.5],'ylim',[0 7],'TickDir','out');
    
    
    subplot(2,2,2);
    ll={}; US=[sigtype{1,:}]; hold on; ii=1;
    for i=1:length(US)
        sel=DT1.signal==categorical(i);
        if sum(sel)>0
            x=[0 MaxECC];
            if i==1
                y=[...
                    mdl_v1.Coefficients.Estimate(1) + ...
                    mdl_v1.Coefficients.Estimate(slidx).*x ...
                    ];
                plot(x,y, 'LineWidth',5, 'Color',CLR(i,:));
                
            else
                j = find(v1idx==i,1,'first');
                y=[...
                    mdl_v1.Coefficients.Estimate(1) + ...
                    mdl_v1.Coefficients.Estimate(1+j) + ...
                    (mdl_v1.Coefficients.Estimate(slidx) + ...
                    mdl_v1.Coefficients.Estimate(slidx+j)).*x ...
                    ];
                plot(x,y, 'LineWidth',3, 'Color',CLR(i,:));
                
            end
            ll{ii}=sigtype{2,i};
            ii=ii+1;
        end
    end
    legend(ll,'Location','NorthEastOutside');
    title('V1 FIT');xlabel('Ecc');ylabel('Sz');
    set(gca,'xlim',[0 12.5],'ylim',[0 7],'TickDir','out');
    
    % v4
    DT4=DT(DT.area==categorical(4) & DT.ecc<=MaxECC,:);
    mdl_v4 = fitlm(DT4,'sz ~ 1 + signal*ecc');
    
    v4idx=[];
    for j=1:length(mdl_v4.CoefficientNames)
        if mdl_v4.CoefficientNames{j}(1:3) == 'sig' & ...
                mdl_v4.CoefficientNames{j}(end-2:end) ~= 'ecc'
            v4idx = [v4idx str2double(mdl_v4.CoefficientNames{j}(8:end))];
        end
    end
    
    subplot(2,2,3);
    ll={}; US=[sigtype{1,:}]; hold on; ii=1;
    for i=1:length(US)
        sel=DT4.signal==categorical(i);
        if sum(sel)>0
            %scatter(DT4.ecc(sel),DT4.sz(sel),'filled','CData',CLR(i,:));
            scatter(DT4.ecc(sel),DT4.sz(sel),msz,...
                'Marker','o',...
                'MarkerEdgeColor',CLR(i,:),'MarkerFaceColor',CLR(i,:),...
                'MarkerEdgeAlpha',0,'MarkerFaceAlpha',0.25);
            ll{ii}=sigtype{2,i};
            ii=ii+1;
        end
    end
    legend(ll,'Location','NorthEastOutside');
    title('V4');xlabel('Ecc');ylabel('Sz');
    set(gca,'xlim',[0 12.5],'ylim',[0 7],'TickDir','out');
    
    subplot(2,2,4);
    ll={}; US=[sigtype{1,:}]; hold on; ii=1;
    for i=1:length(US)
        sel=DT4.signal==categorical(i);
        if sum(sel)>0
            x=[0 MaxECC];
            if i==1
                y=[...
                    mdl_v4.Coefficients.Estimate(1) + ...
                    mdl_v4.Coefficients.Estimate(slidx).*x ...
                    ];
                plot(x,y, 'LineWidth',5, 'Color',CLR(i,:));
            else
                j = find(v4idx==i,1,'first');
                y=[...
                    mdl_v4.Coefficients.Estimate(1) + ...
                    mdl_v4.Coefficients.Estimate(1+j) + ...
                    (mdl_v4.Coefficients.Estimate(slidx) + ...
                    mdl_v4.Coefficients.Estimate(slidx+j)).*x ...
                    ];
                plot(x,y, 'LineWidth',3, 'Color',CLR(i,:));
            end
            ll{ii}=sigtype{2,i};
            ii=ii+1;
        end
    end
    legend(ll,'Location','NorthEastOutside');
    title('V4 FIT');xlabel('Ecc');ylabel('Sz');
    set(gca,'xlim',[0 12.5],'ylim',[0 7],'TickDir','out');
    
    sgtitle(['MODEL: ' MMS{m} ', Rth_mri: ' num2str(Rth_mri) ...
        ', Rth_ephys:' num2str(Rth_ephys)],'interpreter','none')
    if SaveFigs
        saveas(f_xmod,fullfile(figfld, ['XMOD_REGR_' MMS{m} '.png']));
    end
    if CloseFigs; close(f_xmod); end
    
    % Look at effect sizes ===
    %plotEffects(mdl_v1)
    %plotInteraction(mdl_v1,'ecc','signal')
    
    %% WHICH SIGNAL TYPES HAVE SIGNIFICANT ECC-SZ RELATION ====================
    areas = [1 4];
    st=[sigtype{1,:}];sn=sigtype(2,:);
    
    for a = 1:length(areas) % loop over areas
        stats.area(a).signalswithslope = [];
        for s = st % loop over signals
            sel = ...
                DT.area==categorical(areas(a)) & ...
                DT.signal==categorical(s) & ...
                DT.ecc<=MaxECC;
            if sum(sel)>1
                tbl = DT(sel,:);
                mdl = fitlm(tbl,'sz ~ 1 + ecc');
                stats.area(a).signal(s).name = sn{s};
                stats.area(a).signal(s).mdl = mdl;
                stats.area(a).signal(s).anova = anova(mdl);
                stats.area(a).signal(s).ic = mdl.Coefficients(1,:);
                stats.area(a).signal(s).sl = mdl.Coefficients(2,:);
                CI = coefCI(mdl);
                stats.area(a).signal(s).icCI = CI(1,:);
                stats.area(a).signal(s).slCI = CI(2,:);
                if stats.area(a).signal(s).sl.pValue < 0.05
                    stats.area(a).signalswithslope = [...
                        stats.area(a).signalswithslope s];
                end
                
                fprintf(['Area ' num2str(areas(a)) ', Signal: '  num2str(sn{s}) ...
                    ', slope: ' num2str(stats.area(a).signal(s).sl.Estimate) ...
                    ', n = ' num2str(stats.area(a).signal(s).anova.DF(2)+1) ...
                    ', t = ' num2str(stats.area(a).signal(s).sl.tStat) ...
                    ', p = ' num2str(stats.area(a).signal(s).sl.pValue) '\n']);
            end
        end
    end
    
    %% TEST SIGNALS WITH SIGNIFICANT SLOPE TOGETHER IN ONE LINEAR MODEL =======
    for a = 1:length(areas) % loop over areas
        sel = ...
            DT.area==categorical(areas(a)) & ...
            ismember(DT.signal,categorical(stats.area(a).signalswithslope)) & ...
            DT.ecc<=MaxECC;
        
        % recreate a table to avoid empty signal categories
        signal2 = zeros(size(signal_num));
        for i=1:length(stats.area(a).signalswithslope)
            signal2(signal_num==stats.area(a).signalswithslope(i))=i;
        end
        
        signal2 = signal2(sel);
        area2 = area(sel);
        ecc2=ecc(sel);
        sz2=sz(sel);
        
        signal2=categorical(signal2);
        tbl2=table(signal2,area2, ecc2, sz2);
        mdl2 = fitlm(tbl2,'sz2 ~ signal2*ecc2');
        
        stats.area(a).full.mdl = mdl2;
        stats.area(a).full.anova = anova(mdl2);
        stats.area(a).full.signals = sn(stats.area(a).signalswithslope);
        
        fprintf('========================================\n')
        fprintf(['All significant signal types in one model V' num2str(areas(a)) '\n']);
        fprintf('========================================\n')
        fprintf('Eccentricity\n')
        fprintf(['F = ' num2str(stats.area(a).full.anova.F(2)) ...
            ', df = ' num2str(stats.area(a).full.anova.DF(2)) ...
            ', p = ' num2str(stats.area(a).full.anova.pValue(2)) '\n']);
        fprintf('Signal\n')
        fprintf(['F = ' num2str(stats.area(a).full.anova.F(1)) ...
            ', df = ' num2str(stats.area(a).full.anova.DF(1)) ...
            ', p = ' num2str(stats.area(a).full.anova.pValue(1)) '\n']);
        fprintf('Ecc*Signal\n')
        fprintf(['F = ' num2str(stats.area(a).full.anova.F(3)) ...
            ', df = ' num2str(stats.area(a).full.anova.DF(3)) ...
            ', p = ' num2str(stats.area(a).full.anova.pValue(3)) '\n']);
    end
    clear signal2;
    
    %% TEST SIGNALS WITH SIGNIFICANT SLOPE AGAINST MRI ========================
    for a = 1:length(areas) % loop over areas
        ss = stats.area(a).signalswithslope; ss(ss==1)=[];
        sidx=1;
        fprintf('========================================\n')
        fprintf(['MRI vs EPHYS V' num2str(areas(a)) '\n']);
        fprintf('========================================\n')
        for s = ss
            sel = ...
                DT.area == categorical(areas(a)) & ...
                ismember(DT.signal,categorical([1 s])) & ...
                DT.ecc <= MaxECC ;
            
            % recreate a table to avoid empty signal categories
            signal2 = zeros(size(signal_num));
            for i=1:length(stats.area(a).signalswithslope)
                signal2(signal_num==stats.area(a).signalswithslope(i))=i;
            end
            
            signal2 = signal2(sel);
            area2 = area(sel);
            ecc2=ecc(sel);
            sz2=sz(sel);
            
            signal2=categorical(signal2);
            tbl2=table(signal2,area2, ecc2, sz2);
            mdl2 = fitlm(tbl2,'sz2 ~ signal2*ecc2');
            
            stats.area(a).sig_vsMRI(sidx).name = sn{s};
            stats.area(a).sig_vsMRI(sidx).mdl = mdl2;
            stats.area(a).sig_vsMRI(sidx).anova = anova(mdl2);
            
            fprintf(['MRI vs ' stats.area(a).sig_vsMRI(sidx).name '\n']);
            
            fprintf('Eccentricity\n')
            fprintf(['F = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.F(2)) ...
                ', df = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.DF(2)) ...
                ', p = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.pValue(2)) '\n']);
            fprintf('Signal\n')
            fprintf(['F = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.F(1)) ...
                ', df = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.DF(1)) ...
                ', p = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.pValue(1)) '\n']);
            fprintf('Ecc*Signal\n')
            fprintf(['F = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.F(3)) ...
                ', df = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.DF(3)) ...
                ', p = ' num2str(stats.area(a).sig_vsMRI(sidx).anova.pValue(3)) '\n']);
            
            sidx=sidx+1;
        end
    end
    
    %% PLOT SLOPES WITH CI AND INDICATE SIGN DIFF FROM MRI ====================
    f_slopes = figure;
    set(f_slopes,'Position',[100 100 1200 400]);
    
    % V1 -----
    subplot(1,2,1); hold on;
    mri_sl = stats.area(1).signal(1).sl.Estimate;
    mri_sl_bot = stats.area(1).signal(1).slCI(1);
    mri_sl_top = stats.area(1).signal(1).slCI(2);
    x=[0 10]; xarea=[x fliplr(x)];
    
    yarea=[mri_sl_bot mri_sl_bot mri_sl_top mri_sl_top];
    fill(xarea,yarea,'r','EdgeColor','none','FaceAlpha',0.25)
    plot(x,[mri_sl mri_sl],'r','Linewidth',3);
    
    for s=2:7
        NoData=false;
        if isempty(stats.area(1).signal(s).sl)
            NoData=true;
        end
        if~NoData
            Y=[stats.area(1).signal(s).sl.Estimate ...
                stats.area(1).signal(s).slCI(1) ...
                stats.area(1).signal(s).slCI(2)];
            % check if this is significantly different from MRI
            isSim=false;
            for ss = 1:length(stats.area(1).sig_vsMRI)
                if strcmp(stats.area(1).signal(s).name,...
                        stats.area(1).sig_vsMRI(ss).name) & ...
                        stats.area(1).sig_vsMRI(ss).anova.pValue(3) > 0.05
                    isSim=true;
                end
            end
            if isSim
                bar(s,Y(1),'FaceColor',[0.4660, 0.6740, 0.1880]);
                errorbar(s,Y(1),Y(1)-Y(2),Y(3)-Y(1),'k','Linestyle','none');
            elseif NoData
                % nothing to plot
            else
                bar(s,Y(1),'FaceColor',[.25 .25 .25]);
                errorbar(s,Y(1),Y(1)-Y(2),Y(3)-Y(1),'k','Linestyle','none');
            end
        end
    end
    set(gca,'xlim',[1.5 7.5],'ylim',[0 0.5])
    set(gca,'xtick',1:7,'xticklabels',sigtype(2,:))
    ylabel('Ecc-Sz SLOPE');
    title('V1','interpreter','none');
    
    % V4 -----
    subplot(1,2,2); hold on;
    mri_sl = stats.area(2).signal(1).sl.Estimate;
    mri_sl_bot = stats.area(2).signal(1).slCI(1);
    mri_sl_top = stats.area(2).signal(1).slCI(2);
    x=[0 10]; xarea=[x fliplr(x)];
    
    yarea=[mri_sl_bot mri_sl_bot mri_sl_top mri_sl_top];
    fill(xarea,yarea,'r','EdgeColor','none','FaceAlpha',0.25)
    plot(x,[mri_sl mri_sl],'r','Linewidth',3);
    
    for s=2:7
        NoData=false;
        if isempty(stats.area(2).signal(s).sl)
            NoData=true;
        end
        if~NoData
            Y=[stats.area(2).signal(s).sl.Estimate ...
                stats.area(2).signal(s).slCI(1) ...
                stats.area(2).signal(s).slCI(2)];
            % check if this is significantly different from MRI
            isSim=false;
            for ss = 1:length(stats.area(2).sig_vsMRI)
                if strcmp(stats.area(2).signal(s).name,...
                        stats.area(2).sig_vsMRI(ss).name) & ...
                        stats.area(2).sig_vsMRI(ss).anova.pValue(3) > 0.05
                    isSim=true;
                end
            end
            if isSim
                bar(s,Y(1),'FaceColor',[0.4660, 0.6740, 0.1880]);
                errorbar(s,Y(1),Y(1)-Y(2),Y(3)-Y(1),'k','Linestyle','none');
            elseif NoData
                % nothing to plot
            else
                bar(s,Y(1),'FaceColor',[.25 .25 .25]);
                errorbar(s,Y(1),Y(1)-Y(2),Y(3)-Y(1),'k','Linestyle','none');
            end
        end
    end
    set(gca,'xlim',[1.5 7.5],'ylim',[0 0.5])
    set(gca,'xtick',1:7,'xticklabels',sigtype(2,:))
    ylabel('Ecc-Sz SLOPE');
    title('V4','interpreter','none');
    sgtitle(['MODEL: ' MMS{m} ', Rth_mri: ' num2str(Rth_mri) ...
        ', Rth_ephys:' num2str(Rth_ephys)],'interpreter','none')
    if SaveFigs
        saveas(f_slopes,fullfile(figfld, ['XMOD_SlopeComparison_' MMS{m} '.png']));
    end
    if CloseFigs; close(f_slopes); end
end