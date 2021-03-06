addpath(genpath('../../prfCode'));

sess={...
%     '20180201',...
%     '20180131',...
%     '20180125_2',...
%     '20180125_1',...
%     '20180124',...
%     '20180117_2',...
%     '20180117_1',...
%     '20171220',...
%     '20171214_2',...
%     '20171214_1',...
    '20171207',...
    '20171129_2',...
    '20171129_1',...
    '20171116',...
};
for s=1:length(sess)
    pRF_avg_BOLD_tseries('danny',sess{s});
end


sess={...
    '20171129',...
    '20170518',...
    '20170512',...
    '20170411',...
    '20160804',...
    '20160803',...
    '20160729',...
    '20160728',...
    };
for s=1:length(sess)
    pRF_avg_BOLD_tseries('eddy',sess{s});
end

rmpath(genpath('../../prfCode'));