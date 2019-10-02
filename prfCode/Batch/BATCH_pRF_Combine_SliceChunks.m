addpath(genpath('../../prfCode'));

Ms={'danny','eddy'};
Ss={'csshrf_cv1_dhrf',...
    'csshrf_cv1_mhrf',...
    'doghrf_cv1_dhrf','doghrf_cv1_mhrf',...
    'linhrf_cv1_dhrf','linhrf_cv1_mhrf',...
    'linhrf_cv1_mhrf_neggain','linhrf_cv1_dhrf_neggain'...
    };

for m=2%1:length(Ms)
    fprintf(['MONKEY: ' Ms{m} '\n']);
    for s=1:length(Ss)
        pRF_Combine_SliceChunks_cv(Ms{m},Ss{s});
        %pRF_Combine_SliceChunks_cv0(Ms{m},Ss{s});
    end
end

rmpath(genpath('../../prfCode'));