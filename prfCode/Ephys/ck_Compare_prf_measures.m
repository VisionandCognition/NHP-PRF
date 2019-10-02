B1=load('~/Desktop/BOLD_polyfit_danny');
%B2=load('~/Desktop/BOLD_polyfit_eddy');
M=load('~/Desktop/MUA_polyfit');
L=load('~/Desktop/LFP_polyfit');
%%
x=0.5:1:8.5;

figure;
hold on; box on;
y = B1.BOLD_pfit(2)+B1.BOLD_pfit(1)*x;
plot(x,y,'LineWidth',2,'Color',[1 0 0]);
y = B2.BOLD_pfit(2)+B2.BOLD_pfit(1)*x;
plot(x,y,'LineWidth',2,'Color',[1 0.5 0]);
y = M.MUA_pfit(2)+M.MUA_pfit(1)*x;
plot(x,y,'LineWidth',2,'Color','b');
for i=1:5
    y = L.LFP_pfit(i,2)+L.LFP_pfit(i,1)*x;
    plot(x,y,'LineWidth',2,'Color',[0 1/i 0]);
end