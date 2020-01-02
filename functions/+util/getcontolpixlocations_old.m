function [st,ed] = getcontolpixlocations_old(scopeloc, params, scopeparams)
    gridix = scopeloc.gridix ;
    loc = scopeloc.loc ;
    dims = params.imagesize ;
    scopeparams1 = scopeparams(1) ;
    imsize_um = scopeparams1.imsize_um ;
    scopeparams1dims = scopeparams1.dims ;
    targetoverlap_um = 4*[5 5 5]; %in um
    N = params.Ndivs;
    gridix_xyz = gridix(:,1:3) ;
    s1 = median(gridix_xyz);  
    %i1 = find(gridix(:,1)==s1(1)&gridix(:,2)==s1(2)&gridix(:,3)==s1(3)) ;
    i1 = find(all(gridix_xyz == s1, 2)) ;
    s2 = s1+1;
    %i2 = find(gridix(:,1)==s2(1)&gridix(:,2)==s2(2)&gridix(:,3)==s2(3)) ;
    i2 = find(all(gridix_xyz == s2, 2)) ;
    sdiff = abs(diff(loc([i2 i1],:)))*1000
    overlap_um = round(imsize_um-sdiff);
    pixsize = imsize_um./(scopeparams1dims-1);
    ovelap_px = round(overlap_um./pixsize);
    targetoverlap_pix = round(targetoverlap_um./pixsize); %in um
    st = round(ovelap_px/2-targetoverlap_pix/2);
    % find nearest integer that makes ed-st divisible by number of segments (N)
    ed = st+floor((dims-2*st)/N)*N;
end
