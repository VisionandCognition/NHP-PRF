% Compare the linear fits from MS-derived code with that from analyzePRF

%% Load ===================================================================
fprintf('Loading results...\n');
base='/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF';
load(fullfile(base,'FitResults','MultiModal','MUA_Struct'));
load(fullfile(base,'FitResults','MultiModal','old','pRF_estimates_ephys'));
fprintf('Done!\n')

%% set R2 window to look at ===============================================
model='linear_ephys_cv1';
monkey='lick';

RT = [80 100];

%% R2 =====================================================================
figure; 
subplot(2,2,1);hold on;

mR = []; 
for i=1:8
    mR = [mR; RetMap(1).MUA(i).rcrf(1:128)'];
end

aR = MUA.RTEmx.R2(strcmp(...
    MUA.RTEmx.Monkey,monkey) & ...
    strcmp(MUA.RTEmx.Model,model)...
    );

plot([0 100],[0 100],'k');
xlabel('manual fit'); ylabel('analyzePRF');
scatter(mR*100,aR);
axis([0 100 0 100])
title('R2');

%% size ===================================================================
subplot(2,2,2);hold on;

mSz = [];
for i=1:8
    mSz = [mSz; RetMap(1).MUA(i).rsd(1:128)];
end

aSz = MUA.RTEmx.rfs(...
    strcmp(MUA.RTEmx.Monkey,monkey) & ...
    strcmp(MUA.RTEmx.Model,model)...
    );

plot([0 5],[0 5],'k');
xlabel('manual fit'); ylabel('analyzePRF');
scatter(mSz(aR>RT(1) & aR<RT(2)),...
    aSz(aR>RT(1) & aR<RT(2)));
axis([0 5 0 5])
title('size')

dSz = abs(...
    mSz(aR>RT(1) & aR<RT(2)) - ...
    aSz(aR>RT(1) & aR<RT(2))...
    );

%% position ===============================================================
subplot(2,2,3);hold on;

mXY = [];
for i=1:8
    mXY = [mXY; RetMap(1).MUA(i).rmx(1:128) RetMap(1).MUA(i).rmy(1:128)];
end

aX = MUA.RTEmx.X(...
    strcmp(MUA.RTEmx.Monkey,monkey) & ...
    strcmp(MUA.RTEmx.Model,model)...
    );
aY = MUA.RTEmx.Y(...
    strcmp(MUA.RTEmx.Monkey,monkey) & ...
    strcmp(MUA.RTEmx.Model,model)...
    );

plot([-10 10],[-10 10],'k');
xlabel('manual fit'); ylabel('analyzePRF');
scatter(mXY(aR>RT(1) & aR<RT(2),1),...
    aX(aR>RT(1) & aR<RT(2)));
scatter(mXY(aR>RT(1) & aR<RT(2),2),...
    aY(aR>RT(1) & aR<RT(2)));
axis([-10 10 -10 10])
title('X (blue) and Y (red) positions')

dP = sqrt(...
    (...
    mXY(aR>RT(1) & aR<RT(2),1) - ...
    aX(aR>RT(1) & aR<RT(2)) ...
    ).^2 + ...
    (...
    mXY(aR>RT(1) & aR<RT(2),2) - ...
    aY(aR>RT(1) & aR<RT(2)) ...
    ).^2 ...
    );

mECC = sqrt(mXY(:,1).^2 + mXY(:,2).^2 );

aECC = sqrt(aX.^2 + aY.^2);

%% Diff pos vs diff size ==================================================
subplot(2,2,4);hold on;
xlabel('diff pos (abs)'); ylabel('diff sz (abs)');
scatter(dP,dSz);
title('diif size vs diff pos')

%% distributions ==========================================================
figure;
subplot(2,2,1); hold on;
hist(dSz,100); 
title('distibutions of size differences')
subplot(2,2,2); hold on;
hist(dP,100);
title('distibutions of poisiton differences')

%% Size vs ecc ============================================================
%figure;
subplot(2,2,3)
scatter(aECC(aR>RT(1) & aR<RT(2)),aSz(aR>RT(1) & aR<RT(2)))
title('analyzePRF results')
xlabel('Ecc');ylabel('size');
axis([0 10 0 5])
subplot(2,2,4)
scatter(mECC(aR>RT(1) & aR<RT(2)),mSz(aR>RT(1) & aR<RT(2)))
title('manual fits')
xlabel('Ecc');ylabel('size');
axis([0 10 0 5])

%% 
figure; 
subplot(1,2,1);hold on;
s=aR>RT(1) & aR<RT(2);
viscircles([aX(s) aY(s)],aSz(s),'Color','r');
viscircles(mXY(s,:),mSz(s),'Color','b');

subplot(1,2,2);hold on;
for f=find(s==1)'
    plot([mXY(f,1) aX(f)], [mXY(f,2),aY(f)],'o-');
    plot(aX(f), aY(f),'o','MarkerFaceColor','k');
end