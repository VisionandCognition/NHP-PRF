Ms={'danny','eddy'};
Ss={'AllSessions-only_avg','AllSessions-avg-even','AllSessions-avg-odd',...
    'AllSessions-only_avg_DHRF','AllSessions-avg-even_DHRF','AllSessions-avg-odd_DHRF'};
for m=1:length(Ms)
    fprintf(['MONKEY: ' Ms{m} '\n']);
    for s=1:length(Ss)
        pRF_Combine_SliceChunks(Ms{m},Ss{s});
    end
end