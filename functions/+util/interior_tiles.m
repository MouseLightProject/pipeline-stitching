function interior = interior_tiles(scopeloc,num_dilate)
%%
if nargin<2
    num_dilate=0;
end
grids = scopeloc.gridix(:,1:3);
grids = grids-min(grids)+1;
dims = range(grids)+1;
dims = dims([2 1 3]);
tileImage = zeros(dims);

inds = sub2ind(dims,grids(:,2),grids(:,1),grids(:,3));
tileImage(inds) = 1;
out = bwperim(tileImage);
for it = 1:num_dilate
    out = imdilate(out,strel(ones(3,3,3))) .* tileImage;
end

out_inds = find(out);
M = containers.Map(inds,1:size(grids,1));

these_inds = zeros(size(out_inds));
for ii=1:size(out_inds,1)
    these_inds(ii) = M(out_inds(ii));
end

interior = zeros(size(grids,1),1);
interior(these_inds) = 1;
interior = ~interior;
end
