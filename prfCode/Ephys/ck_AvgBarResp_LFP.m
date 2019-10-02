function RPB = ck_AvgBarResp_LFP(a,c,L,EVT,pos)
%a=array
%c=channel
%L=L response

%% sweep starts
sw_start=[]; rst_start=[];
for p=1:length(pos)
    if EVT.StimName{pos(p)}==0
        rst_start=[rst_start;pos(p)];
        np = EVT.t_corr(pos(p));
    elseif mod(EVT.StimName{pos(p)},30)==1
        sw_start=[sw_start;pos(p)];
        np = EVT.t_corr(pos(p));
    end
end
RPB.sw_start=sw_start;
RPB.rst_start=rst_start;

%% average over runs
runstart=[];runtraces=[];
for p=1:length(pos)
    if EVT.StimName{pos(p)}==1
        np = EVT.t_corr(pos(p));
        runstart=[runstart np];
    end
end
nsamp=ceil((runstart(2)-runstart(1))*500); % LFP @ 500 Hz

% start 2s before stim onset
runtraces=[];
for r=1:length(runstart)
    si = find(L(a).t >= runstart(r)-2,1,'first');
    runtraces=[runtraces;L(a).chan{c}.data(si:si+nsamp)'];
    if r==1
        mt=L(a).t(si:si+nsamp);
    end
end
mResp=mean(runtraces,1);
BL=mean(mResp(mt<runstart(1)));
RPB.mt=mt;
RPB.mResp=mResp;
RPB.BL=BL;


% plot the stim moments
bi=[];avg_barresp=[];mBarResp=[];
for p=1:length(pos)
    np = EVT.t_corr(pos(p));
    if np<max(mt)
        pp=find(mt>=np,1,'first');
        bi=[bi; pp pp+500];
        mbr=mean(mResp(pp:pp+500));
        mBarResp=[mBarResp;mbr];
        avg_barresp=[avg_barresp; ...
           mt(pp:pp+500)' mbr*ones(501,1)];
    end
end
RPB.mBarResp=mBarResp;
RPB.avg_barresp=avg_barresp;
RPB.bi=bi;

 


