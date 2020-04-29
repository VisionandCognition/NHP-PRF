function ck_AvgEphysSessions(MONKEY,SESS,TYPE)

%clear all; 
close all; clc;

if nargin < 2
    MONKEY = {'Lick','Aston'};
    SESS={'20180807_B2','20191203_B1','20191203_B2','20191203_B3';...
          '20181004_B1','20191205_B1','20191205_B2','20191205_B3'}; % all
    %SESS={'20191203_B2','20191203_B3','20191203_B3';...
    %      '20191205_B1','20191205_B2','20191205_B3'}; % only the new ones
    TYPE='ALL'; % AVG / NEW /ALL
end

data_fld = '/media/NETDISKS/VS02/VandC/PRF_EPHYS/Data_proc';

%%
for m=1:length(MONKEY)
    fprintf(['- ' MONKEY{m} ,' ---\n']);
    res_fld = fullfile(data_fld, MONKEY{m},TYPE);
    mua_res_fld = fullfile(res_fld,'MUA');
    lfp_res_fld = fullfile(res_fld,'LFP'); 
    
    [~,~]=mkdir(mua_res_fld);
    [~,~]=mkdir(lfp_res_fld);

    %% MUA ================================================================
    fprintf('Averaging MUA\n');
    for I=1:8 % instances
        fprintf(['Instance ' num2str(I) '\n']);
        for elec=1:128
            col_bar{elec}=[]; col_bar_even{elec}=[]; col_bar_odd{elec}=[]; 
        end
        
        for sess = 1:size(SESS,2)
            load(fullfile(data_fld,MONKEY{m},SESS{m,sess},'MUA',...
                [MONKEY{m} '_' SESS{m,sess} '_array_' num2str(I) '_mMUA' ]))
            load(fullfile(data_fld,MONKEY{m},SESS{m,sess},'MUA',...
                [MONKEY{m} '_' SESS{m,sess} '_array_' num2str(I) '_mMUA_even' ]),'mMUA_even')
            load(fullfile(data_fld,MONKEY{m},SESS{m,sess},'MUA',...
                [MONKEY{m} '_' SESS{m,sess} '_array_' num2str(I) '_mMUA_odd' ]),'mMUA_odd')
            
            for elec=1:128
                col_bar{elec}=[col_bar{elec}; mMUA(elec).bar]; 
                col_bar_even{elec}=[col_bar_even{elec}; mMUA_even(elec).bar];
                col_bar_odd{elec}=[col_bar_odd{elec}; mMUA_odd(elec).bar];
            end
            clear 'mMUA' 'mMUA_even' 'mMUA_odd' 
        end
        
        if strcmp(TYPE,ALL)
            save(fullfile(mua_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_TRACES']),...
                'col_bar','col_bar_even','col_bar_odd')
        else
            % average
            fprintf('Saving...');
            for elec=1:128
                mMUA(elec).bar=mean(col_bar{elec},1);
            end
            fprintf('ALL...');
            save(fullfile(mua_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_mMUA']),'mMUA','C')
            
            for elec=1:128
                mMUA_even(elec).bar=mean(col_bar_even{elec},1);
            end
            fprintf('EVEN...');
            save(fullfile(mua_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_mMUA_even']),'mMUA_even','C')
            
            for elec=1:128
                mMUA_odd(elec).bar=mean(col_bar_odd{elec},1);
            end
            fprintf('ODD...\n');
            save(fullfile(mua_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_mMUA_odd']),'mMUA_odd','C')
        end
    end
    clear 'col_bar' 'col_bar_even' 'col_bar_odd'
    
    %% LFP ================================================================
    fprintf('Averaging LFP\n');
    for I=1:8 % instances
        fprintf(['Instance ' num2str(I) '\n']);
        for elec=1:128
            for fb = 1:5
                col_bar{elec,fb}=[]; 
                col_bar_even{elec,fb}=[]; 
                col_bar_odd{elec,fb}=[];
                
                col_bl{elec,fb}=[]; 
                col_bl_even{elec,fb}=[]; 
                col_bl_odd{elec,fb}=[];
            end
        end
        
        for sess = 1:size(SESS,2)
            load(fullfile(data_fld,MONKEY{m},SESS{m,sess},'LFP',...
                [MONKEY{m} '_' SESS{m,sess} '_array_' num2str(I) '_mLFP' ]))
            load(fullfile(data_fld,MONKEY{m},SESS{m,sess},'LFP',...
                [MONKEY{m} '_' SESS{m,sess} '_array_' num2str(I) '_mLFP_even' ]),'mLFP_even')
            load(fullfile(data_fld,MONKEY{m},SESS{m,sess},'LFP',...
                [MONKEY{m} '_' SESS{m,sess} '_array_' num2str(I) '_mLFP_odd' ]),'mLFP_odd')
            
            for elec=1:128
                for fb = 1:5
                    col_bar{elec,fb}=[col_bar{elec,fb}; ...
                        mLFP(elec).freq(fb).bar - mLFP(elec).freq(fb).BL]; 
                    col_bar_even{elec,fb}=[col_bar_even{elec,fb}; ...
                        mLFP_even(elec).freq(fb).bar - mLFP_even(elec).freq(fb).BL];
                    col_bar_odd{elec,fb}=[col_bar_odd{elec,fb}; ...
                        mLFP_odd(elec).freq(fb).bar - mLFP_odd(elec).freq(fb).BL];
                    
                    col_bl{elec,fb}=[col_bl{elec,fb}; 0]; 
                    col_bl_even{elec,fb}=[col_bl_even{elec,fb}; 0];
                    col_bl_odd{elec,fb}=[col_bl_odd{elec,fb}; 0];
                end
            end
            clear 'mLFP' 'mLFP_even' 'mLFP_odd' 
        end
        
        if strcmp(TYPE,ALL)
            save(fullfile(lfp_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_TRACES']),...
                'col_bar','col_bar_even','col_bar_odd')
        else
            % average
            fprintf('Saving...');
            for elec=1:128
                for fb = 1:5
                    mLFP(elec).freq(fb).bar=mean(col_bar{elec,fb},1);
                    mLFP(elec).freq(fb).BL=mean(col_bl{elec,fb},1);
                end
            end
            fprintf('ALL...');
            save(fullfile(lfp_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_mLFP']),'mLFP','C')
            
            for elec=1:128
                for fb = 1:5
                    mLFP_even(elec).freq(fb).bar=mean(col_bar_even{elec,fb},1);
                    mLFP_even(elec).freq(fb).BL=mean(col_bl_even{elec,fb},1);
                end
            end
            fprintf('EVEN...');
            save(fullfile(lfp_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_mLFP_even']),'mLFP_even','C')
            
            for elec=1:128
                for fb = 1:5
                    mLFP_odd(elec).freq(fb).bar=mean(col_bar_odd{elec,fb},1);
                    mLFP_odd(elec).freq(fb).BL=mean(col_bl_odd{elec,fb},1);
                end
            end
            fprintf('ODD...\n');
            save(fullfile(lfp_res_fld,...
                [MONKEY{m} '_AVG_array_' num2str(I) '_mLFP_odd']),'mLFP_odd','C')
        end
    end
    clear 'col_bar' 'col_bar_even' 'col_bar_odd'
end