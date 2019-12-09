% Compare the linear fits from MS-derived code with that from analyzePRF

%% Load ===================================================================
CVMODE='cv1';
fprintf('Loading results...\n');
base='/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF';
load(fullfile(base,'FitResults','MultiModal',CVMODE,'LFP_Struct'));
load(fullfile(base,'FitResults','MultiModal','old','pRF_estimates_ephys'));
fprintf('Done!\n')

%% set R2 window to look at ===============================================
model=['linear_ephys_' CVMODE];
monkey='lick';

RT = [50 100];
maxSz = 1000;
SigTypes = {'Theta','Alpha','Beta','lGamma','hGamma'};

%%
for  fb=1:5
    %% R2 =================================================================
    mR = []; 
    for i=1:8
        mR = [mR; RetMap(1).LFP(i,fb).rcrf(1:128)'];
    end

    aR = LFP.RTEmx.R2(...
        strcmp(LFP.RTEmx.Monkey,monkey) & ...
        strcmp(LFP.RTEmx.Model,model) & ...
        strcmp(LFP.RTEmx.SigType,SigTypes{fb}) ...
        );
    
    figure;
    subplot(1,2,1);
    histogram(mR*100,0:100);
    title(['manual (' SigTypes{fb} ')']);
    subplot(1,2,2);
    histogram(aR,0:100);
    title(['analyzePRF (' SigTypes{fb} ')']);
    
    figure;
    subplot(2,2,1);hold on;
    
    plot([0 100],[0 100],'k');
    xlabel('manual fit'); ylabel('analyzePRF');
    scatter(mR*100,aR);
    axis([0 100 0 100])
    title(['R2 (' SigTypes{fb} ')']);
    
    %% size ===============================================================
    subplot(2,2,2);hold on;
    
    mSz = [];
    for i=1:8
        mSz = [mSz; RetMap(1).LFP(i,fb).rsd(1:128)];
    end

    aSz = LFP.RTEmx.rfs(...
        strcmp(LFP.RTEmx.Monkey,monkey) & ...
        strcmp(LFP.RTEmx.Model,model) & ...
        strcmp(LFP.RTEmx.SigType,SigTypes{fb}) ...
        );

    plot([0 5],[0 5],'k');
    xlabel('manual fit'); ylabel('analyzePRF');
    scatter(mSz(aR>RT(1) & aR<RT(2)),...
        aSz(aR>RT(1) & aR<RT(2)));
    axis([0 5 0 5])
    title(['Size (' SigTypes{fb} ')'])

    dSz = abs(...
        mSz(aR>RT(1) & aR<RT(2)  & aSz<maxSz) - ...
        aSz(aR>RT(1) & aR<RT(2)  & aSz<maxSz)...
        );

    %% position ===========================================================
    subplot(2,2,3);hold on;

    mXY = [];
    for i=1:8
        mXY = [mXY; RetMap(1).LFP(i,fb).rmx(1:128) ...
            RetMap(1).LFP(i,fb).rmy(1:128)];
    end
    
    aX = LFP.RTEmx.X(...
        strcmp(LFP.RTEmx.Monkey,monkey) & ...
        strcmp(LFP.RTEmx.Model,model) & ...
        strcmp(LFP.RTEmx.SigType,SigTypes{fb}) ...
        );
    aY = LFP.RTEmx.Y(...
        strcmp(LFP.RTEmx.Monkey,monkey) & ...
        strcmp(LFP.RTEmx.Model,model) & ...
        strcmp(LFP.RTEmx.SigType,SigTypes{fb}) ...
        );
    
    plot([-10 10],[-10 10],'k');
    xlabel('manual fit'); ylabel('analyzePRF');
    scatter(mXY(aR>RT(1) & aR<RT(2) & aSz<maxSz,1),...
        aX(aR>RT(1) & aR<RT(2) & aSz<maxSz));
    scatter(mXY(aR>RT(1) & aR<RT(2) & aSz<maxSz,2),...
        aY(aR>RT(1) & aR<RT(2) & aSz<maxSz));
    axis([-10 10 -10 10])
    title(['X (blue) and Y (red) positions (' SigTypes{fb} ')'])
    
    dP = sqrt(...
        (...
        mXY(aR>RT(1) & aR<RT(2) & aSz<maxSz,1) - ...
        aX(aR>RT(1) & aR<RT(2) & aSz<maxSz) ...
        ).^2 + ...
        (...
        mXY(aR>RT(1) & aR<RT(2) & aSz<maxSz,2) - ...
        aY(aR>RT(1) & aR<RT(2) & aSz<maxSz) ...
        ).^2 ...
        );

    mECC = sqrt(mXY(:,1).^2 + mXY(:,2).^2 );
    aECC = sqrt(aX.^2 + aY.^2);

    %% Diff pos vs diff size ==============================================
    subplot(2,2,4);hold on;
    xlabel('diff pos (abs)'); ylabel('diff sz (abs)');
    scatter(dP,dSz);
    title(['Diff size vs Diff pos (' SigTypes{fb} ')'])

    %% distributions ======================================================
    figure;
    subplot(2,2,1); hold on;
    hist(dSz,100);
    title(['Distributions of size diff (' SigTypes{fb} ')'])
    subplot(2,2,2); hold on;
    hist(dP,100);
    title(['Distributions of position diff (' SigTypes{fb} ')'])

    %% Size vs ecc ========================================================
    subplot(2,2,3)
    scatter(aECC(aR>RT(1) & aR<RT(2)),aSz(aR>RT(1) & aR<RT(2)))
    title(['analyzePRF results (' SigTypes{fb} ')'])
    xlabel('Ecc');ylabel('size');
    axis([0 10 0 5])
    subplot(2,2,4)
    scatter(mECC(aR>RT(1) & aR<RT(2)),mSz(aR>RT(1) & aR<RT(2)))
    title(['manual fits (' SigTypes{fb} ')'])
    xlabel('Ecc');ylabel('size');
    axis([0 10 0 5])

%% 
    figure; 
    subplot(1,3,1);hold on;
    s=aR>RT(1) & aR<RT(2);
    set(gca,'xlim',[-5 15],'ylim',[-15 5]);
    viscircles([aX(s) aY(s)],aSz(s),'Color','r');
    title(['analyzeprf (' SigTypes{fb} ')'])
    subplot(1,3,2);hold on;
    set(gca,'xlim',[-5 15],'ylim',[-15 5]);
    viscircles(mXY(s,:),mSz(s),'Color','b');
    title(['manual (' SigTypes{fb} ')'])

    subplot(1,3,3);hold on;
    for f=find(s==1)'
        plot([mXY(f,1) aX(f)], [mXY(f,2),aY(f)],'o-');
        plot(aX(f), aY(f),'o','MarkerFaceColor','k');
    end
    title(['LFP (' SigTypes{fb} ')'])
end