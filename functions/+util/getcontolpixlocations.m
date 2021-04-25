function [st,ed] = getcontolpixlocations(scopeloc, params, scopeparams)
    gridix = scopeloc.gridix ;
    loc = scopeloc.loc ;
    dims = params.imagesize ;
    scopeparams1 = scopeparams(1) ;
    imsize_um = scopeparams1.imsize_um ;
    scopeparams1dims = scopeparams1.dims ;
    targetoverlap_um = 4*[5 5 5]; %in um
    N = params.Ndivs;
    gridix_xyz = gridix(:,1:3) ;
    %s1 = median(gridix_xyz);  
    %i1 = find(gridix(:,1)==s1(1)&gridix(:,2)==s1(2)&gridix(:,3)==s1(3)) ;
    %i1 = find(all(gridix_xyz == s1, 2)) ;
    %s2 = s1+1;
    %i2 = find(gridix(:,1)==s2(1)&gridix(:,2)==s2(2)&gridix(:,3)==s2(3)) ;
    %i2 = find(all(gridix_xyz == s2, 2)) ;
    %sdiff_old = abs(diff(loc([i2 i1],:)))*1000 ;
    sdiff = abs(compute_sdiff(gridix_xyz, loc))*1000 ;
    overlap_um = round(imsize_um-sdiff);
    pixsize = imsize_um./(scopeparams1dims-1);
    ovelap_px = round(overlap_um./pixsize);
    targetoverlap_pix = round(targetoverlap_um./pixsize); %in um
    st = round(ovelap_px/2-targetoverlap_pix/2);
    % find nearest integer that makes ed-st divisible by number of segments (N)
    ed = st+floor((dims-2*st)/N)*N;
end



function result = compute_sdiff(gridix_xyz, loc)
    % Robust determination of the spacing between tiles
    dimension_count = size(gridix_xyz, 2) ;
    result = nan(1, dimension_count) ;
    for i = 1 : dimension_count ,
        result(i) = compute_sdiff_1d(gridix_xyz(:,i), loc(:,i)) ;
    end
end



function result = compute_sdiff_1d(is, xs)
    levels = unique(is) ;
    level_count = length(levels) ;
    xs_from_level_index = nan(size(levels)) ;
    for level_index = 1 : level_count ,
        level = levels(level_index) ;
        is_this_level = (is == level) ;
        if any(is_this_level) ,
            xs_for_this_level = xs(is_this_level) ;
            xs_from_level_index(level_index) = median(xs_for_this_level) ;
        end
    end    
    interlevel_deltas = diff(xs_from_level_index) ;
    result = nanmedian(interlevel_deltas) ;  %#ok<NANMEDIAN> % median, ignoring nans
end
