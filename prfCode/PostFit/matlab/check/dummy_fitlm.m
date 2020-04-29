x=1:100;
% 1
IC=0; sl=1;
y1=IC + sl*x + 5*(rand(1,100)-0.5); % slope only (intercept=0)
% 2
IC=4; sl=1;
y2=IC + sl*x + 5*(rand(1,100,1)-0.5); % other intercept, same slope
% 3 
IC=0; sl=1.5;
y3=IC + sl*x + 5*(rand(1,100)-0.5); % same intercept, other slope (small diff)
% 4
IC=0; sl=5;
y4=IC + sl*x + 5*(rand(1,100)-0.5); % same intercept, other slope (larger diff)
% 5
IC=5; sl=8;
y5=IC + sl*x + 5*(rand(1,100)-0.5); % other intercept, other slope


%
g=categorical([ones(1,100) 2*ones(1,100) 3*ones(1,100) 4*ones(1,100) 5*ones(1,100)]');
x=[x x x x x]'; y=[y1 y2 y3 y4 y5]';

tbl=table(x,y,g);

xx=x(101:200); yy=y2'; gg=g(101:200);
tbl2=table(xx,yy,gg);

mdl = fitlm(tbl,'y ~ 1 + g*x')
mdl2 = fitlm(tbl2,'yy ~ 1 + gg*xx')

% cov matrix
% i1   sl1    i2    sl2    i3    sl3    i4    sl4

ct=zeros(1,10); ct(1,2)=0; ct(1,4)=1;
[p,f,r] = coefTest(mdl,ct)

anova(mdl)