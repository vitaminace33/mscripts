function [V,DV,GV] = harmosc(q,order)
narginchk(1,2)
if nargin == 1
    order = 0;
end
nargoutchk(0,3-order)
nout = max(1,nargout);

switch order
    case 0
        V = q.^2/2;
        if nout > 1
            DV = q;
            if nout > 2
                GV = q;
            end
        end
    case 1
        V = q;
        if nout > 1
            DV = q;
        end
    case 2
        V = q;
end
end
