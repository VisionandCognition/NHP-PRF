function out = ck_2dPRF_ephys(xc,yc,sd, settings)

% creates a 3d array of 2d pRF heatmaps
% xc : x-coordinate(s) in degrees (Nx1 vector)
% yc : y-coordinate(s) in degrees (Nx1 vector)
% sd : sima(s) in degrees (Nx1 vector)
% settings.PixPerDeg : pixels per dva
% settings.meshsize : 1000

if nargin < 4
    settings.PixPerDeg = 29.5032;
    settings.meshsize = 2000;   
else
    if ~isfield(settings,'PixPerDeg')
        settings.PixPerDeg = 29.5032;
    end
    if ~isfield(settings,'meshsize')
        settings.meshsize = 2000;
    end
end

xmesh = ((1:settings.meshsize)-settings.meshsize/2)./settings.PixPerDeg;
ymesh = ((settings.meshsize:-1:1)-settings.meshsize/2)./settings.PixPerDeg;

% [x y sd]
gdets = [xc,yc,sd];  
ngauss = size(gdets,1);

ZZ=[];
for gz = 1:ngauss
    %Read in Gaussian parameters
    a = gdets(gz,1);
    b = gdets(gz,2);
    c = gdets(gz,3);
    %Make Gaussian
    gw = normpdf(xmesh,a,c);
    gv = normpdf(ymesh,b,c);
    Z = gv'*gw;
    Z = Z./sum(sum(Z)); %Could normalise by sum
    if gz == 1
        ZZ = zeros(size(Z,1),size(Z,2),ngauss,'single');
    end
    ZZ(:,:,gz) = single(Z);
end

out.img = ZZ;
out.xmesh = xmesh;
out.ymesh = ymesh;
out.xr = [xmesh(1) xmesh(end)];
out.yr = [ymesh(end) ymesh(1)];
