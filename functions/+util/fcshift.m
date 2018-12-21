function [locs,xshift2D,yshift2D] = fcshift(model,order,xy,dims,locs)
% estimates 2D FC fields as a function of curvature model and distance to
% imaging center. 
if isempty(xy)
    xlocs = 1:dims(1);
    ylocs = 1:dims(2);
    [xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
    xy = [xy1(:),xy2(:)];
end
cent = squeeze(mean(model(1:2,1),3));
scale = (squeeze(mean(model(1:2,2),3)));
shift = squeeze(mean(model(1:2,3),3));
beta = scale./shift.^order;
[xshift2D,yshift2D] = util.shiftxy(xy,cent,beta,order,dims);
idxctrl = sub2ind(dims([2 1]),locs(:,2),locs(:,1));
xshift = xshift2D(idxctrl);
yshift = yshift2D(idxctrl);
locs(:,1) = locs(:,1) + xshift;
locs(:,2) = locs(:,2) + yshift;
end