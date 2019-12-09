%% look at full timecourse relative to stimuli
a=1; % array
c=1; % channel
f1=figure; hold on;
% trace
sw = 100;% in ms org data is 1kHz
plot(L(a).t(L(a).t>0),smooth(L(a).chan{c}.data(L(a).t>0),sw));

%% stim starts
ypos = median(L(a).chan{c}.data(L(a).t>0));
for p=1:length(pos)
    np = EVT.t_corr(pos(p));
    plot([np np],[ypos-2 ypos+2],'r')
end

%% sweep starts
sw_start=[]; rst_start=[];
for p=1:length(pos)
    if EVT.StimName{pos(p)}==0
        rst_start=[rst_start;pos(p)];
        np = EVT.t_corr(pos(p));
        plot([np np],[ypos-2 ypos+2],'g','LineWidth',2)
    elseif mod(EVT.StimName{pos(p)},30)==1
        sw_start=[sw_start;pos(p)];
        plot([np np],[ypos-2 ypos+2],'r')
        np = EVT.t_corr(pos(p));
        plot([np np],[ypos-2 ypos+2],'k','LineWidth',2)
    end
end
close(f1)

%% average over runs
runstart=[];
runtraces=[];
for p=1:length(pos)
    if EVT.StimName{pos(p)}==1
        np = EVT.t_corr(pos(p));
        runstart=[runstart np];
    end
end
nsamp=ceil((runstart(2)-runstart(1))*500);

% start 2s before stim onset
for r=1:length(runstart)
    si = find(L(a).t >= runstart(r)-2,1,'first');
    runtraces=[runtraces;L(a).chan{c}.data(si:si+nsamp)'];
    if r==1
        mt=L(a).t(si:si+nsamp);
    end
end
mResp=mean(runtraces,1);
BL=mean(mResp(mt<runstart(1)));

% plot the averaged trace
f2=figure; hold on;
sw = 100;% in ms org data is 1kHz
plot(mt,smooth(mResp,sw)-BL);

% plot the stim moments
bi=[];avg_barresp=[];mBarResp=[];
for p=1:length(pos)
    np = EVT.t_corr(pos(p));
    if np<max(mt)
        plot([np np],[-1 2],'r')
        pp=find(mt>=np,1,'first');
        bi=[bi; pp pp+500];
        mbr=mean(mResp(pp:pp+500));
        mBarResp=[mBarResp; mt(pp+250) mbr];
        avg_barresp=[avg_barresp; ...
           mt(pp:pp+500)' mbr*ones(501,1)];
    end
end

% plot the sweep starts
for p=1:length(pos)
    if mod(EVT.StimName{pos(p)},30)==1
        np = EVT.t_corr(pos(p));
        if np<max(mt)
            plot([np np],[-1 2],'k','LineWidth',2)
        end
    end
end
close(f2)

%% digitize to a single response per bar
% plot the averaged trace
f3=figure; hold on;
bar(mBarResp(:,1),mBarResp(:,2)-BL,'BarWidth',1)
maxY=1.2*max(mBarResp(:,2)-BL);
minY=1.2*min(mBarResp(:,2)-BL);

% plot the stim moments
for p=1:length(pos)
    np = EVT.t_corr(pos(p));
    if np<max(mt)
        plot([np np],[minY maxY],'r')
        pp=find(mt>=np,1,'first');
    end
end

% plot the sweep starts
for p=1:length(pos)
    if mod(EVT.StimName{pos(p)},30)==1
        np = EVT.t_corr(pos(p));
        if np<max(mt)
            plot([np np],[minY maxY],'k','LineWidth',2)
        end
    end
end
 


