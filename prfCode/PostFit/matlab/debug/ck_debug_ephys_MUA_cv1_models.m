% Compare the linear fits from MS-derived code with that from analyzePRF

%% Load ===================================================================
CVMODE={'cv1','cv1'};
DS={'ORG','NEW'};
fprintf('Loading results...\n');
base='/Users/chris/Dropbox/CURRENT_PROJECTS/NHP_MRI/Projects/pRF';
load(fullfile(base,'FitResults','MultiModal',DS{1},CVMODE{1},'MUA_Struct'));
MUA1=MUA;
load(fullfile(base,'FitResults','MultiModal',DS{2},CVMODE{2},'MUA_Struct'));
MUA2=MUA;
fprintf('Done!\n')
xvalsel = 'RTEmx';

%% set R2 window to look at ===============================================
model1=['linear_ephys_' CVMODE{1}];
model2=['linear_ephys_' CVMODE{2}];

monkey='aston';

RT = [0 100];
maxSz = 1000;

%% R2 =====================================================================
aR1 = MUA1.(xvalsel).R2(strcmp(...
    MUA1.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA1.(xvalsel).Model,model1)...
    );

aR2 = MUA2.(xvalsel).R2(strcmp(...
    MUA2.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA2.(xvalsel).Model,model2)...
    );

figure;
subplot(1,2,1);
histogram(aR1,0:100);
title(CVMODE{1});
subplot(1,2,2);
histogram(aR2,0:100);
title(CVMODE{2});

figure; 
subplot(2,2,1);hold on;

plot([0 100],[0 100],'k');
xlabel(CVMODE{1}); ylabel(CVMODE{2});
scatter(aR1,aR2);
axis([0 100 0 100])
title('R2');

%% size ===================================================================
subplot(2,2,2);hold on;

aSz1 = MUA1.(xvalsel).rfs(...
    strcmp(MUA1.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA1.(xvalsel).Model,model1)...
    );

aSz2 = MUA2.(xvalsel).rfs(...
    strcmp(MUA2.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA2.(xvalsel).Model,model2)...
    );

%% Universal selection ====================================================
S0 = (aR1>RT(1) & aR1<RT(2) & aR2>RT(1) & aR2<RT(2));
S = (aR1>RT(1) & aR1<RT(2) & aR2>RT(1) & aR2<RT(2) & aSz1<maxSz);
LRF = aSz1<maxSz;

dSz = abs( aSz1(S) - aSz2(S) );

plot([0 5],[0 5],'k');
xlabel(CVMODE{1}); ylabel(CVMODE{2});
scatter(aSz1(S0),aSz2(S0));
axis([0 5 0 5])
title('size')

%% position ===============================================================
subplot(2,2,3);hold on;

aX1 = MUA1.(xvalsel).X(...
    strcmp(MUA1.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA1.(xvalsel).Model,model1)...
    );
aY1 = MUA1.(xvalsel).Y(...
    strcmp(MUA1.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA1.(xvalsel).Model,model1)...
    );

aX2 = MUA2.(xvalsel).X(...
    strcmp(MUA2.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA2.(xvalsel).Model,model2)...
    );
aY2 = MUA2.(xvalsel).Y(...
    strcmp(MUA2.(xvalsel).Monkey,monkey) & ...
    strcmp(MUA2.(xvalsel).Model,model2)...
    );

plot([-10 10],[-10 10],'k');
xlabel(CVMODE{1}); ylabel(CVMODE{2});
scatter(aX1(S),aX2(S));
scatter(aY1(S),aY2(S));
axis([-10 10 -10 10])
title('X (blue) and Y (red) positions')

dP = sqrt((aX1(S) - aX2(S) ).^2 + (aY1(S) - aY2(S) ).^2 );

aECC1 = sqrt(aX1.^2 + aY1.^2);
aECC2 = sqrt(aX2.^2 + aY2.^2);

%% Diff pos vs diff size ==================================================
subplot(2,2,4);hold on;
xlabel('diff pos (abs)'); ylabel('diff sz (abs)');
scatter(dP,dSz);
title('diif size vs diff pos')

%% distributions ==========================================================
figure;
subplot(2,2,1); hold on;
hist(dSz,100); 
title('distributions of size differences')
subplot(2,2,2); hold on;
hist(dP,100);
title('distributions of position differences')

%% Size vs ecc ============================================================
%figure;
subplot(2,2,3)
scatter(aECC1(S0),aSz1(S0))
title(CVMODE{1})
xlabel('Ecc');ylabel('size');
axis([0 10 0 5])
subplot(2,2,4)
scatter(aECC2(S0),aSz2(S0))
title(CVMODE{2})
xlabel('Ecc');ylabel('size');
axis([0 10 0 5])

%% Position vs eccentricity ===============================================
figure;
subplot(2,2,1);
scatter(aX1(S0),aSz1(S0))
title(CVMODE{1})
xlabel('X');ylabel('size');
axis([-10 10 0 35])
subplot(2,2,2)
scatter(aX2(S0),aSz2(S0))
title(CVMODE{2})
xlabel('X');ylabel('size');
axis([-10 10 0 35])

subplot(2,2,3);
scatter(aY1(S0),aSz1(S0))
title(CVMODE{1})
xlabel('Y');ylabel('size');
axis([-15 5 0 35])
subplot(2,2,4)
scatter(aY2(S0),aSz2(S0))
title(CVMODE{2})
xlabel('Y');ylabel('size');
axis([-15 5 0 35])

%% PRF in VF ==============================================================
figure; 
subplot(1,2,1);hold on;
s=aR1>RT(1) & aR1<RT(2);
set(gca,'xlim',[-5 15],'ylim',[-15 5]);
title(CVMODE{1})
viscircles([aX1(s) aY1(s)],aSz1(s),'Color','r');
subplot(1,2,2);hold on;
viscircles([aX2(s) aY2(s)],aSz2(s),'Color','b');
set(gca,'xlim',[-5 15],'ylim',[-15 5]);
title(CVMODE{2})

%% Distance between different model based pRFs in VF ======================
figure;hold on;
for f=find(s==1)'
    plot([aX1(f) aX2(f)], [aY1(f),aY2(f)],'o-');
    plot(aX2(f), aY2(f),'o','MarkerFaceColor','k');
end