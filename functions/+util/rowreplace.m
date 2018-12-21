function X = rowreplace(X,inds,row)
% replaces certain rows of X based on inds array with the value in row.
% faster form of: for ix = inds(:)';X(ix,:) = row;end
for i = 1:size(X,2)
    X(inds,i) = row(i);
end