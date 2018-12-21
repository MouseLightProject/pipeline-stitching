function [st,ed] = getcontolpixlocations(scopeloc,params,scopeparams)

dims = params.imagesize;
targetoverlap_um = 4*[5 5 5]; %in um
N = params.Ndivs;
s1 = median(scopeloc.gridix(:,1:3));  i1=find(scopeloc.gridix(:,1)==s1(1)&scopeloc.gridix(:,2)==s1(2)&scopeloc.gridix(:,3)==s1(3));
s1 = s1+1;i2=find(scopeloc.gridix(:,1)==s1(1)&scopeloc.gridix(:,2)==s1(2)&scopeloc.gridix(:,3)==s1(3));
sdiff = abs(diff(scopeloc.loc([i2 i1],:)))*1000;
overlap_um = round(scopeparams(1).imsize_um-sdiff);
pixsize = scopeparams(1).imsize_um./[scopeparams(1).dims-1];
ovelap_px = round(overlap_um./pixsize);
targetoverlap_pix = round(targetoverlap_um./pixsize); %in um
st = round(ovelap_px/2-targetoverlap_pix/2);
% find nearest integer that makes ed-st divisible by number of segments (N)
ed = st+floor((dims-2*st)/N)*N;

