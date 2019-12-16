function check_newephys

%clear all; 
clc

%% load the ephys data
M={'Lick','Aston'};
fld='ALL';

%% MUA ====================================================================
fprintf('Creating plots of MUA traces...');
signal='MUA';
for m=1:length(M)
    data_fld = fullfile('/media/NETDISKS/VS02/VandC/PRF_EPHYS/Data_proc',...
        M{m},fld,signal);
    save_fld = ['EphysCheck_' M{m} '_' fld '_' signal];
    [~,~]=mkdir(save_fld);

    for I=1:8 % instances
        load(fullfile(data_fld,[M{m} '_AVG_array_' num2str(I) '_TRACES'])) %#ok<*LOAD>
        for elec=1:128
            spi=mod(elec,16);
            if spi==0; spi=16; end

            if spi==1
                f=figure('visible','off');hold on;
                set(f,'Position',[0 0 1800 1000]);
                first_elec = elec;
            end
            
            subplot(4,4,spi)
            plot(col_bar{elec}','Linewidth',1) %#ok<*IDISVAR,*USENS>
            set(gca,'xticklabels',{},'yticklabels',{},...
                'xaxislocation','origin')
            legend({'ORG','New1','New2','New3'},...
                'Orientation','Horizontal','Location','BestOutside')
            title(['Inst ' num2str(I) ', Ch ' num2str(elec)])
        
            if spi==16
                saveas(f,fullfile(save_fld,...
                    ['Inst-' num2str(I) '_Elec-' ...
                    num2str(first_elec,'%03.f') 'to' num2str(elec,'%03.f') '.png']));
                close(f);
            end
        end
    end
end
fprintf('DONE\n');

%% LFP ====================================================================
fprintf('Creating plots of LFP traces...');
signal='LFP'; fb='hGAM';
for m=1:length(M)
    data_fld = fullfile('/media/NETDISKS/VS02/VandC/PRF_EPHYS/Data_proc',...
        M{m},fld,signal);
    save_fld = ['EphysCheck_' M{m} '_' fld '_' signal];
    [~,~]=mkdir(save_fld);

    for I=1:8 % instances
        load(fullfile(data_fld,[M{m} '_AVG_array_' num2str(I) '_TRACES']))
        for elec=1:128
            spi=mod(elec,16);
            if spi==0; spi=16; end

            if spi==1
                f=figure('visible','off');hold on;
                set(f,'Position',[0 0 1800 1000]);
                first_elec = elec;
            end
            
            subplot(4,4,spi)
            plot(col_bar{elec,5}','Linewidth',1)
            set(gca,'xticklabels',{},'yticklabels',{},...
                'xaxislocation','origin')
            legend({'ORG','New1','New2','New3'},...
                'Orientation','Horizontal','Location','BestOutside')
            title(['Inst ' num2str(I) ', Ch ' num2str(elec)])
        
            if spi==16
                saveas(f,fullfile(save_fld,...
                    ['Inst-' num2str(I) '_Elec-' ...
                    num2str(first_elec,'%03.f') 'to' num2str(elec,'%03.f') ...
                    '_' fb '.png']));
                close(f);
            end
        end
    end
end
fprintf('DONE\n');
