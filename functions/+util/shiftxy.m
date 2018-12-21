function [xshift,yshift] = shiftxy(xy,centxy,beta,order,dims)
if nargin<4
    order = 1;
end
% we split contribution to half as curvature is symetric
% beta: p(2)/p(3), weight: x-xcent
% at edge of a tile beta*weight will be roughly p(2)/2 as p(3) is stage
% shift and x-xcent will will be around dims/2. for order =1, this results
% in warp at edge half of the estimated curvature:
% beta*weight = p(2)/dim_image * (dim_image-dim_image/2) = p(2)/2
% beta*weight*(x-p(1))^2 = p(2)*(x-p(1))^2 / 2 = curvature_model/2

repthis = (isnan(centxy) | centxy==0);
centxy(repthis) = dims(repthis);

weightx = ((xy(:,1)-centxy(2)).^order);
weighty = ((xy(:,2)-centxy(1)).^order);

xshift = beta(1)*weightx.*((xy(:,2)-centxy(1)).^2);
yshift = beta(2)*weighty.*((xy(:,1)-centxy(2)).^2);

if nargin>4
    xshift = reshape(xshift,dims([2 1]));
    yshift = reshape(yshift,dims([2 1]));
end
