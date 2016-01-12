function [mu, sigma, n_ind, l, r, gamma, Gamma] = mces(x,varargin)

if nargin > 1
    seq = varargin(1);
else
    seq = 'monotone';
end

mu = mean(x);
x = x - mu;

n = length(x);
m = floor(n/2);

gamma = zeros(1,n);
Gamma = zeros(1,m);

for i = 1:m
    j1 = 2*i;
    j0 = j1-1;
    gamma(j0) = sum(x(1:n-j0+1).*x(j0:n))/n;
    gamma(j1) = sum(x(1:n-j0  ).*x(j1:n))/n;
    Gamma(i) = gamma(j0) + gamma(j1);
    if Gamma(i) <= 0
        break;
    end
end

m = i-1;
Gamma = Gamma(1:m);

if any(strcmp(seq,{'pos','positive'}))
    % done
elseif any(strcmp(seq,{'mono','monotone'}))
    Gamma = gdm(Gamma);
elseif any(strcmp(seq,{'conv','convex'}))
    Gamma = gcm(gdm(Gamma));
end

sigma2 = gamma(1);
sigma  = sqrt(sigma2);
Sigma2 = -sigma2+2*sum(Gamma);
if sigma2 > 0
    n_ind  = round(n/(Sigma2/sigma2)); % Sigma2/n = sigma2/n_ind
else
    n_ind = 1;
end

coef = 1.96; % si 95%
%coef = 2.576; % si 99%
shift = coef*sqrt(Sigma2/n);
l = mu-shift;
r = mu+shift;
end

% Greatest deacreasing minorant
function x = gdm(x)
for i = 2:length(x)
    x(i) = min(x(i-1),x(i));
end
end

% Greatest convex minorant
function x = gcm(x)
m = length(x);
i = 1;
while i < m
    [slope,j] = min([(x(i+1:m)-x(i))./(1:m-i),-x(i)/(m+1-i)]);
    x(i+1:i+j-1) = (1:j-1)*slope + x(i);
    i = i + j;
end
end
