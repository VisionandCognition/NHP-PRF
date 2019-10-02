clear all; clc;
add_classicRF = true; % don't have to do this every time. Once suffices.

MONKEY={'Lick','Aston'};
data_fld=pwd;
inst=1:8;

load('pRF_estimates_ephys')

if add_classicRF
    for m=1:length(MONKEY)
        for a=inst
            fprintf(['Loading ' MONKEY{m} ', inst' num2str(a) '\n']);
            M(m).inst(a) = load(fullfile(data_fld,...
                ['classicRF' MONKEY{m}],['RFs_instance' num2str(a)]),...
                'RFs','channelRFs','meanChannelSNR');
        end
    end
    
    Pix2Deg = M(1).inst(1).RFs{1}.szdeg / M(1).inst(1).RFs{1}.sz;
    
    for m=1:length(MONKEY)
        table_arr = [];
        VarNames = {...
            'Instance', 'ChanIndex', 'Array', ...
            'Chan', 'Area', ...
            'rf_x_pix', 'rf_x_deg', ...
            'rf_y_pix','rf_y_deg', ...
            'rf_sz_pix', 'rf_sz_deg', ...
            'rf_ang', 'rf_theta', ...
            'rf_ecc', 'rf_ecc_deg', ...
            'rf_hrad_pix', 'rf_hrad_deg', ...
            'rf_vrad_pix', 'rf_vrad_deg', ...
            'mSNR', ...
            };
        for a=inst
            RetMap(m).classic_RF = M(m).inst(a).RFs;
            RetMap(m).classic_RF_snr = M(m).inst(a).meanChannelSNR;
            RetMap(m).classic_RF_chan = M(m).inst(a).channelRFs;
            
            for c=1:length(M(m).inst(a).RFs)
                table_arr = [table_arr;...
                    a c RetMap(m).ChanMap.arrayNums(c,a) ...
                    RetMap(m).ChanMap.channelNums(c,a) RetMap(m).ChanMap.areas(c,a) ...
                    M(m).inst(a).RFs{c}.centrex M(m).inst(a).RFs{c}.centrex*Pix2Deg ...
                    M(m).inst(a).RFs{c}.centrey M(m).inst(a).RFs{c}.centrey*Pix2Deg ...
                    M(m).inst(a).RFs{c}.sz M(m).inst(a).RFs{c}.szdeg ...
                    M(m).inst(a).RFs{c}.ang M(m).inst(a).RFs{c}.theta ...
                    M(m).inst(a).RFs{c}.ecc M(m).inst(a).RFs{c}.ecc*Pix2Deg ...
                    M(m).inst(a).RFs{c}.horRad M(m).inst(a).RFs{c}.horRad*Pix2Deg ...
                    M(m).inst(a).RFs{c}.verRad M(m).inst(a).RFs{c}.verRad*Pix2Deg ...
                    M(m).inst(a).meanChannelSNR(c) ...
                    ];
            end
        end
        RetMap(m).classic_RF_table = array2table(table_arr,'VariableNames', VarNames);
    end
    fprintf('Resaving pRF_estimates_ephys.mat with classic RF definition added\n');
    save('pRF_estimates_ephys','RetMap')
end

%% plot some things
for v = [1 4]
    f1=figure; f2=figure;
    for m=1:2
        figure(f1);
        areas = RetMap(1).table_mua.Area;
        
        rf_x = RetMap(m).classic_RF_table.rf_x_deg;
        prf_x = RetMap(m).table_mua.prf_x;
        rf_y = RetMap(m).classic_RF_table.rf_y_deg;
        prf_y = RetMap(m).table_mua.prf_y;
        rf_sz = RetMap(m).classic_RF_table.rf_sz_deg;
        prf_sz = RetMap(m).table_mua.prf_sd*2;
        
        rf_ecc = RetMap(m).classic_RF_table.rf_ecc_deg;
        prf_ecc= RetMap(m).table_mua.prf_ecc;
        
        rf_snr = RetMap(m).classic_RF_table.mSNR;
        prf_r2 = RetMap(m).table_mua.prf_r2;
        
        TH = [1 0.50]; % [min_SNR min_R2]
        
        idx = rf_snr>=TH(1) & prf_r2>=TH(2) & areas == v;
        
        subplot(3,2,0+m); hold on;
        plot([-10 10],[-10 10],'-r')
        plot(rf_x(idx),prf_x(idx),'ok')
        title(['V' num2str(v) ' ' MONKEY{m} ' x'])
        set(gca,'xlim',[-2 10],'ylim',[-2 10])
        xlabel('RF');ylabel('PRF');
        [r,p] = corr(rf_x(idx),prf_x(idx));
        text(-1,9,num2str(r))
        text(-1,8,num2str(p))
        
        subplot(3,2,2+m); hold on;
        plot([-10 10],[-10 10],'-r')
        plot(rf_y(idx),prf_y(idx),'ok')
        title(['V' num2str(v) ' ' MONKEY{m} ' y'])
        set(gca,'xlim',[-10 5],'ylim',[-10 5])
        xlabel('RF');ylabel('PRF');
        [r,p] = corr(rf_y(idx),prf_y(idx));
        text(-9,4,num2str(r))
        text(-9,3,num2str(p))
        
        subplot(3,2,4+m); hold on;
        plot([-10 10],[-10 10],'-r')
        plot(rf_sz(idx),prf_sz(idx),'ok')
        title(['V' num2str(v) ' ' MONKEY{m} ' sz'])
        set(gca,'xlim',[0 10],'ylim',[0 10])
        xlabel('RF');ylabel('PRF');
        [r,p] = corr(rf_sz(idx),prf_sz(idx));
        text(1,9,num2str(r))
        text(1,8,num2str(p))
        
        figure(f2);
        subplot(1,2,m); hold on
        plot(rf_ecc(idx),rf_sz(idx),'ok')
        plot(prf_ecc(idx),prf_sz(idx),'or')
        title(['V' num2str(v) ' ' MONKEY{m} ' ecc vs sz'])
        set(gca,'xlim',[0 10],'ylim',[0 10])
        xlabel('ecc');ylabel('sz');
    end
end

%% look per array
colind = hsv(16);
colind(8,:)=[139/255 69/255 19/255];
% select area
v=[1 4]; % if vector all areas mentioned will be included
f1=figure;
for m=1:2
    figure(f1);
    areas = RetMap(1).table_mua.Area;
    arr =  RetMap(1).table_mua.Array;
    
    rf_x = RetMap(m).classic_RF_table.rf_x_deg;
    prf_x = RetMap(m).table_mua.prf_x;
    rf_y = RetMap(m).classic_RF_table.rf_y_deg;
    prf_y = RetMap(m).table_mua.prf_y;
    rf_sz = RetMap(m).classic_RF_table.rf_sz_deg;
    prf_sz = RetMap(m).table_mua.prf_sd*2;
    
    rf_ecc = RetMap(m).classic_RF_table.rf_ecc_deg;
    prf_ecc= RetMap(m).table_mua.prf_ecc;
    
    rf_snr = RetMap(m).classic_RF_table.mSNR;
    prf_r2 = RetMap(m).table_mua.prf_r2;
    
    TH = [1 0.70]; % [min_SNR min_R2]
    LT = cell(1,16);
    LT{1}='X=Y';
    for a=1:16
        idx = rf_snr>=TH(1) & prf_r2>=TH(2) & arr == a & ismember(areas, v);
        
        subplot(3,2,0+m); hold on;
        if a==1; plot([-10 10],[-10 10],'-r'); end
        plot(rf_x(idx),prf_x(idx),'o','MarkerEdgeColor',colind(a,:))
        title([MONKEY{m} ' x'])
        set(gca,'xlim',[-2 10],'ylim',[-2 10])
        xlabel('RF');ylabel('PRF');
        
        subplot(3,2,2+m); hold on;
        if a==1; plot([-10 10],[-10 10],'-r'); end
        plot(rf_y(idx),prf_y(idx),'o','MarkerEdgeColor',colind(a,:))
        title([MONKEY{m} ' y'])
        set(gca,'xlim',[-10 5],'ylim',[-10 5])
        xlabel('RF');ylabel('PRF');
        
        
        subplot(3,2,4+m); hold on;
        if a==1; plot([-10 10],[-10 10],'-r'); end
        plot(rf_sz(idx),prf_sz(idx),'o','MarkerEdgeColor',colind(a,:))
        title([MONKEY{m} ' sz'])
        set(gca,'xlim',[0 10],'ylim',[0 10])
        xlabel('RF');ylabel('PRF');
        LT{a+1}=['Arr ' num2str(a)];
    end
end
subplot(3,2,3);
l=legend(LT); set(l,'Location','NorthWest');
legend