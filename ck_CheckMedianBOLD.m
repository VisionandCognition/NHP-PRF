imagesc(medianBOLD(:,:,44,1))

%%
v1=[56 15 44];
v2=[38 16 44];

B1=squeeze([medianBOLD(v1(1),v1(2),v1(3),:)]);
%B2=squeeze([medianBOLD(v2(1),v2(2),v2(3),:)]);
%B1=squeeze([sessmedianBOLD(v1(1),v1(2),v1(3),:)]);
%B2=squeeze([sessmedianBOLD(v2(1),v2(2),v2(3),:)]);

%figure;
hold on;
plot(B1);
plot(B1-median(B1));
%plot(B2-median(B2));
set(gca,'xlim',[0 410]);



%%



imagesc(sessmedianBOLD(:,:,24,1))

%%

B1=squeeze([sessmedianBOLD(v1(1),v1(2),v1(3),:)]);
B2=squeeze([sessmedianBOLD(v2(1),v2(2),v2(3),:)]);

figure;
hold on;
plot(B1-median(B1));
plot(B2-median(B2));
