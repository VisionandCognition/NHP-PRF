Ms={'danny','eddy'};
Ss={'AllSessions-only_avg','AllSessions-avg-even','AllSessions-avg-odd'};
for m=1:length(Ms)
    fprintf(['MONKEY: ' Ms{m} '\n']);
    for s=1:length(Ss)
        pRF_Combine_SliceChunks(Ms{m},Ss{s});
    end
end