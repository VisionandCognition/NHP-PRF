% compressive spatial summation model (Kay & Winnawer)=====================
% 5 parameters

% center gauss
x  = pp(1)
y  = pp(2) 
sd = pp(3)
gain = pp(4) % (amplitude is set to 1, this gain will be applied applied to both Gs)

% redundant
expt = pp(5) % >>> not used in DOG

% new parameters
sdgain   = pp(5) % sd2/sd1 constrained to be >1
norm_amp = pp(6) % a2/a1 where a is set to always be 1  (later scaled by gain)

% surr gauss
x  = pp(1) % same for both Gs
y  = pp(2) % same for both Gs
sd = pp(3)/pp(5) % sdratio*sd_G1
amp = pp(6) % norm_amp/1


% old linear model
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
    makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
    (2*pi*abs(pp(3))^2))); 0]),options.hrf,dd(:,prod(res)+1);

%conv2run(a,b,c)
% <a> is a 2D matrix with time x cases
posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),...
    makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ...
    (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5))
% posrect() >> x<0 = 0  , pp(4) is gain parameter
posrect(pp(4))  % GAIN

% gaussian 1 (center, positive)
makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2) % gaussian 1
%>>> Here we need to create something new for gaussian 2 (surround,
%negative)
-  pp(6) .* makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)/pp(5)),abs(pp(3)/pp(5)),xx,yy,0,0) / (2*pi*abs(pp(3)/pp(5))^2) % gaussian 2


% so the new DOG model becomes
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res), ...
    (makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2)) - ...
    (pp(6) .* makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)/pp(5)),abs(pp(3)/pp(5)),xx,yy,0,0) / (2*pi*abs(pp(3)/pp(5))^2)) ...
    )) ; 0]),options.hrf,dd(:,prod(res)+1));



%% TO DO
- create supergridseeds for 6 parameters when doing DOG
- 