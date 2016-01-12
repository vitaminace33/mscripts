function [m,n] = linfit(x,y)
  sz = size(y);
  if length(x)==numel(x) && length(x)==sz(1)
      x = repmat(x(:),[1 sz(2:end)]);
  end
  rsz = @(A) repmat(A,[sz(1) 1]);
  m = dot(x,y-rsz(mean(y)))./dot(x,x-rsz(mean(x)));
  n = mean(y-rsz(m).*x);
end
