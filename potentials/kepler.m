function [V,DV,GV] = kepler(q,order)
narginchk(1,2)
if nargin == 1
    order = 0;
end
nargoutchk(0,3-order)
nout = max(1,nargout);

d_2 = 1/(q(1)^2+q(2)^2);
d_1 = sqrt(d_2);
  V = -d_1;
 DV =  d_1*d_2*q;
if nout + order == 3
    d_6 = d_2^3;
     GV = -2*d_6*q;
end

switch order
    case 0
        % do nothing
    case 1
        V = DV;
        if nout == 2
            DV = GV;
        end
    case 2
        V = GV;
end
end
