function [locs_out,xshift2D,yshift2D] = fcshift(model,order,xy,dims,locs_in)
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
    [xshift2D, yshift2D] = shiftxy(xy, cent, scale, shift, order, dims) ;
    idxctrl = sub2ind(dims([2 1]),locs_in(:,2),locs_in(:,1));
    xshift = xshift2D(idxctrl);
    yshift = yshift2D(idxctrl);
    dimension_count = size(locs_in,2) ;
    if dimension_count==2 ,
        locs_shift = [xshift yshift] ;
    elseif dimension_count==3 ,
        location_count = size(locs_in,1) ;
        zshift = zeros(location_count,1) ;
        locs_shift = [xshift yshift zshift] ;
    else
        error('Unsupported dimension count (%d)', dimension_count)
    end
    locs_out = locs_in + locs_shift ;
end



function [xshift, yshift] = shiftxy(xy, cent, scale, shift, order, dims)
    % we split contribution to half as curvature is symetric
    % beta: p(2)/p(3), weight: x-xcent
    % at edge of a tile beta*weight will be roughly p(2)/2 as p(3) is stage
    % shift and x-xcent will will be around dims/2. for order =1, this results
    % in warp at edge half of the estimated curvature:
    % beta*weight = p(2)/dim_image * (dim_image-dim_image/2) = p(2)/2
    % beta*weight*(x-p(1))^2 = p(2)*(x-p(1))^2 / 2 = curvature_model/2
    beta = scale./shift.^order;
    
    repthis = (isnan(cent) | cent==0);
    cent(repthis) = dims(repthis);
    
    weightx = ((xy(:,1)-cent(2)).^order);
    weighty = ((xy(:,2)-cent(1)).^order);
    
    xshift = beta(1)*weightx.*((xy(:,2)-cent(1)).^2);
    yshift = beta(2)*weighty.*((xy(:,1)-cent(2)).^2);
    
    if nargin>4
        xshift = reshape(xshift,dims([2 1]));
        yshift = reshape(yshift,dims([2 1]));
    end
end
