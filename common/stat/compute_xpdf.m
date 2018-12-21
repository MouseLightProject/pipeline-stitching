function [x2, n, b] = compute_xpdf(x)
  x2 = reshape(x, 1, numel(x));
  [n, b] = histcounts(x2);
  % This is definitely not probability density function
  x2 = sort(x2);
  % downsampling to speed up computations
  if numel(x) > 1e4
      x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));
  elseif numel(x) > 10
      x2 = interp1 (1:length(x2), x2, 1:length(x2)/10:length(x2));
  else
      x2 = interp1 (1:length(x2), x2, 1:length(x2));
  end
end