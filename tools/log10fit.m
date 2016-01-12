function [m,n] = log10fit(x,y)
  [m,n] = linfit(log10(x),log10(y));
end
