function RPB = ck_AvgBarResp(a,c,M,EVT,pos)
%a=array
%c=channel
%M=mua response

%% sweep starts
sw_start=[]; rst_start=[];
for p=1:length(pos)
    if EVT.StimName{pos(p)}==0
        rst_start=[rst_start;pos(p)];
        np = EVT.t_corr(pos(p));
        %plot([np np],[ypos-2 ypos+2],'g','LineWidth',2)
    elseif mod(EVT.StimName{pos(p)},30)==1
        sw_start=[sw_start;pos(p)];
        %plot([np np],[ypos-2 ypos+2],'r')
        np = EVT.t_corr(pos(p));
        %plot([np np],[ypos-2 ypos+2],'k','LineWidth',2)
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
nsamp=ceil((runstart(2)-runstart(1))*1000);

% start 2s before stim onset
for r=1:length(runstart)
    si = find(M(a).tsec >= runstart(r)-2,1,'first');
    runtraces=[runtraces;M(a).chan{c}(si:si+nsamp)'];
    if r==1
        mt=M(a).tsec(si:si+nsamp);
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
        %plot([np np],[-1 2],'r')
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

 


