function [toppts,xg,yg] = getCtrlPts(xpix,ypix, params,st,ed)
% calculate the locations, in pixels, of the control points (same for every
% tile).  This depends on the tile dimensions and the number of subvolumes
% to create
N = params.Ndivs;
if (xpix/N)~=round(xpix/N) || (ypix/N)~=round(ypix/N)
    error(['Image dimensions not divisible by ' num2str(N)]);
end

if nargin<4
    xlen = xpix/N;
    ylen = ypix/N;
    xg = 0:xlen-1:xpix;
    yg = 0:ylen-1:ypix;
else
    sdif = ed-st;
    xlen = sdif(1)/(N);
    ylen = sdif(2)/(N);
    xg = st(1):xlen:ed(1);
    yg = st(2):ylen:ed(2);
end

[xgrid, ygrid] = meshgrid(xg, yg);
toppts(:,1) = xgrid(:);
toppts(:,2) = ygrid(:);
%toppts(:,3) = 0;

[~, IX] = sort(toppts(:,2), 'ascend');
toppts = toppts(IX,:);
end