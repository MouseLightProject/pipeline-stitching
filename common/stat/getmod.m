function mod = getmod(x)
% assumes x[Nxdim]
if ndims(x) ==1
    mod = gtmd(x);
else
    if size(x,1)<=size(x,2)
        error('input has more columns then rows, use getmod(x'') if input is correct')
    end
    mod = size(x,2);
    for dim = 1:size(x,2)
        mod(dim) = gtmd(x(:,dim));
    end
end
end

function md = gtmd(x)
    [n,edges] = histcounts(x);
    [maxn,indmax] = max(n);
    md = (edges(indmax)+edges(indmax+1))/2;
end