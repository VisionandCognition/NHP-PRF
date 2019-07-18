Ms={'danny','eddy'};
Ss={'AllSessions-avg-cv','AllSessions-avg-cv_DHRF'};
for m=1:length(Ms)
    fprintf(['MONKEY: ' Ms{m} '\n']);
    for s=1:length(Ss)
        pRF_Combine_SliceChunks_cv(Ms{m},Ss{s});
    end
end