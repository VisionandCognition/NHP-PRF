addpath(genpath('../../FitResults/Ephys'));

load('pRF_estimates_all')

rth=0.7;
AreaName = {0,1,4;'All','V1','V4'};

for m = 1:2
    for VA = 1:3
        % visual field coverage
        f(m,1)=figure;hold on;
        f(m,2)=figure;
        for j=1:2
            subplot(1,2,j);hold on
            plot([0 0],[-100 100],'k','linewidth',.5);
            plot([-100 100],[0 0],'k','linewidth',.5);
        end
        
        for a=[1:8]
            if AreaName{1,VA} % not zero
                inc = RetMap(m).table_mua.Instance == a & RetMap(m).table_mua.prf_r2 > rth & ...
                    RetMap(m).table_mua.Area == AreaName{1,VA};
            else 
                inc = RetMap(m).table_mua.Instance == a & RetMap(m).table_mua.prf_r2 > rth;
            end
            AN = AreaName{2,VA};
            
            figure(f(m,1));
            subplot(2,4,a); hold on; box on;
            plot([0 0],[-100 100],'k','linewidth',.5);
            plot([-100 100],[0 0],'k','linewidth',.5);
            viscircles([RetMap(m).table_mua.prf_x(inc),RetMap(m).table_mua.prf_y(inc)],...
                RetMap(m).table_mua.prf_sd(inc)./2,'Color','r','LineWidth',1);
            text(11.5,-10.5,['R2 > ' num2str(rth)],...
                'HorizontalALignment','right','FontSize',14)
            text(11.5,2.5,['n = ' num2str(sum(inc))],...
                'HorizontalALignment','right','FontSize',14)
            
            set(gca,'xlim',[-4 12],'ylim',[-12 4],'FontSize',14);
            title([RetMap(m).Subj ': Instance ' num2str(a) ', Area: ' AN ' [MUA]'],'FontSize', 18);
            
            figure(f(m,2));
            subplot(1,2,1);hold on;
            plot(RetMap(m).table_mua.prf_x(inc),RetMap(m).table_mua.prf_y(inc),...
                '*','MarkerSize',8);
            set(gca,'xlim',[-4 12],'ylim',[-12 4],'FontSize',14);
            title([RetMap(m).Subj ': All, Area: ' AN ' [MUA]'],'FontSize', 18);
            
            subplot(1,2,2);hold on;
            viscircles([RetMap(m).table_mua.prf_x(inc),RetMap(m).table_mua.prf_y(inc)],...
                RetMap(m).table_mua.prf_sd(inc)./2,'LineWidth',1);
            set(gca,'xlim',[-4 12],'ylim',[-12 4],'FontSize',14);
            title([RetMap(m).Subj ': All, Area: ' AN ' [MUA]'],'FontSize', 18);
        end
        set(f(m,1),'Position',[0 0 1200 600])
        set(f(m,2),'Position',[0 0 1200 600])
        figure(f(m,2));
        for j=1:2
            subplot(1,2,j);
            if AreaName{1,VA} % not zero
                inc = RetMap(m).table_mua.prf_r2 > rth & ...
                    RetMap(m).table_mua.Area == AreaName{1,VA};
            else
                inc = RetMap(m).table_mua.prf_r2 > rth;
            end
            text(11.5,-10.5,['R2 > ' num2str(rth)],...
                'HorizontalALignment','right','FontSize',14)
            text(11.5,2.5,['n = ' num2str(sum(inc))],...
                'HorizontalALignment','right','FontSize',14)
        end
    end
end

rmpath(genpath('../../FitResults/Ephys'));
