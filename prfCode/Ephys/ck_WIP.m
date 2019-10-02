rth=.8
figure;
for m=1:2
    subplot(1,2,m); hold on;
    v1_rfsz{m} = RetMap(m).table_mua.prf_sd(RetMap(m).table_mua.Area==1 & ...
        RetMap(m).table_mua.prf_r2>rth & RetMap(m).table_mua.prf_sd>0.2);
    v4_rfsz{m} = RetMap(m).table_mua.prf_sd(RetMap(m).table_mua.Area==4 & ...
        RetMap(m).table_mua.prf_r2>rth & RetMap(m).table_mua.prf_sd>0.2);
    bins=0:0.5:10;
    B1 = hist(v1_rfsz{m},bins);
    B4 = hist(v4_rfsz{m},bins);
    bar(bins-0.25,B1./max(B1),'BarWidth',1);
    bar(bins-0.25,B4./max(B4),'BarWidth',1);
    
    plot( [mean(v1_rfsz{m}) mean(v1_rfsz{m})],[0 1.25]);
    plot( [mean(v4_rfsz{m}) mean(v4_rfsz{m})],[0 1.25]);
end

%%
figure;
for m=1:2
    subplot(1,2,m); hold on;
    v1_rfsz{m} = RetMap(m).table_mua.prf_sd(RetMap(m).table_mua.Area==1 & RetMap(m).table_mua.prf_r2>0.7);
    v1_rfecc{m} = RetMap(m).table_mua.prf_ecc(RetMap(m).table_mua.Area==1 & RetMap(m).table_mua.prf_r2>0.7);
    
    v4_rfsz{m} = RetMap(m).table_mua.prf_sd(RetMap(m).table_mua.Area==4 & RetMap(m).table_mua.prf_r2>0.7);
    v4_rfecc{m} = RetMap(m).table_mua.prf_ecc(RetMap(m).table_mua.Area==4 & RetMap(m).table_mua.prf_r2>0.7);
    
    e1=[];s1=[];
    e4=[];s4=[];
    for ee=0.5:0.5:15
        idx=v1_rfecc{m}<ee & v1_rfecc{m}>ee-0.5;
        if sum(idx)>0
            e1=[e1 ee-0.25];
            s1=[s1 mean(v1_rfsz{m}(idx))];
        end
        idx=v4_rfecc{m}<ee & v4_rfecc{m}>ee-0.5;
        if sum(idx)>0
            e4=[e4 ee-0.25];
            s4=[s4 mean(v4_rfsz{m}(idx))];
        end
    end
    %[p1,s1]=polyfit(v1_rfecc{m},v1_rfsz{m},1);
    %[p4,s4]=polyfit(v4_rfecc{m},v4_rfsz{m},1);
    
    scatter(v1_rfecc{m},v1_rfsz{m});
    scatter(v4_rfecc{m},v4_rfsz{m});
    
    %plot([0 20], p1(2)+p1(1)*([0 20]))
    %plot([0 20], p4(2)+p4(1)*([0 20]))
    plot(e1,s1,'o');
    plot(e4,s4,'o');
end

idx = find(...
    RetMap(2).table_mua.Area==1 & ...
    RetMap(2).table_mua.prf_ecc>3 & ...
    RetMap(2).table_mua.prf_sd<1 & ...
    RetMap(2).table_mua.prf_r2>0.7);

RetMap(2).table_mua.Array(idx)

