addpath(genpath('../../prfCode'));

Ms={'danny','eddy'};

Ss={'csshrf_cv1_dhrf','csshrf_cv1_mhrf',...
    'doghrf_cv1_dhrf','doghrf_cv1_mhrf',...
    'linhrf_cv1_dhrf','linhrf_cv1_mhrf',...
    'linhrf_cv1_dhrf_neggain','linhrf_cv1_mhrf_neggain'...
    };

for m=1:length(Ms)
    fprintf(['MONKEY: ' Ms{m} '\n']);
    for s=1:length(Ss)
        ck_Combine_MRI_SliceChunks_cv(Ms{m},Ss{s});
    end
end

rmpath(genpath('../../prfCode'));