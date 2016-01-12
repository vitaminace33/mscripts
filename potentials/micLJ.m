function [V,DV,HVv] = micLJ(q,params,order)
% MICLJ computes de Lennard-Jones potential (and derivative tensors) of a
% set of particles using the minimum image convention.
%
%   [1] M.P. Allen, D.J. Tildesley: "Computer Simulation of Liquids",
%       Oxford Science Publications (1989).
%
%   URL: http://google.com/+cedricmartinezcampos
%   URL: http://github.com/vitaminace33
%
%   Copyleft 2014-2015 Cedric M. Campos, Dr., joint work with
%                      J.M. Sanz-Serna, Prof. Dr.
%
%   This file is part of my "mscripts" GitHub repository.
%
%   "mscripts" is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   Foobar is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
%
%   Supported by the European Union under the European Social Fund and the
%   Junta de Castilla y Leon (Spain).
%
%   $Revision: 1.0.0.00 $  $Date: 2016/01/11 18:37:51 $

narginchk(2,3)
if nargin == 2
    order = 0;
end
nargoutchk(0,3-order)
nout = max(1,nargout);

epsilon = params.epsilon;
  sigma = params.sigma;
 sigma6 = sigma^6;
      L = params.lattice;
     L2 = L/2;
  dcut2 = (params.cut)^2;
if isfield(params,'mass')
     mass = params.mass;
end

[dim, n] = size(q);
if dim ~= 3
    if n == 1
        iscol = true;
        n = dim/3;
        dim = 3;
        q = reshape(q,[dim n]);
    else
        error('integrators:micLJ:wrongsize', ...
              'The input q must be of size [3 mols] or [3*mols 1].')
    end
else
    iscol = false;
end
N = n*(n-1)/2;

switch nout + order

    case 1

 V = 0;
for i = 1:n-1
    r = mod(q(:,i+1:n)-q(:,i)*ones(1,n-i)+L2,L)-L2;
    for j = i+1:n
        ji = j - i;
%         d2 = 0;                   %
%         for k = 1:dim             % d2 = norm(r(:,ji))^2;
%             d2 = d2 + r(k,ji)^2;  % for loop faster (<) than norm
%         end                       % for < norm < '* << dot
        d2 = r(1,ji)^2 + r(2,ji)^2 + r(3,ji)^2;
        if d2 < dcut2
            d4 = d2^2;
            sd6 = sigma6/(d4*d2);
            V = V + sd6*(sd6-1);
        end
    end
end
V =  4*epsilon*V;

    case 2

 V = 0;
DV = zeros([dim n]);
for i = 1:n-1
    r = mod(q(:,i+1:n)-q(:,i)*ones(1,n-i)+L2,L)-L2;
    for j = i+1:n
        ji = j - i;
%         d2 = 0;                   %
%         for k = 1:dim             % d2 = norm(r(:,ji))^2;
%             d2 = d2 + r(k,ji)^2;  % for loop faster (<) than norm
%         end                       % for < norm < '* << dot
        d2 = r(1,ji)^2 + r(2,ji)^2 + r(3,ji)^2;
        if d2 < dcut2
            d4 = d2^2;
            sd6 = sigma6/(d4*d2);
            V = V + sd6*(sd6-1);
            Dv = sd6*(1-2*sd6)/d2;
            Dvr = Dv*r(:,ji);
            DV(:,i) = DV(:,i) - Dvr;
            DV(:,j) = DV(:,j) + Dvr;
        end
    end
end
 V =  4*epsilon*V;
DV = 24*epsilon*DV;

if iscol
    DV = DV(:);
end

    case 3

persistent index %#ok<TLEV>             % hack: lowers computation, however
if isempty(index) || size(index,2)~=n;  %   increases slightly memory load.
    index = zeros([n-1 n]);
    for i = 1:n-1
        for j = i+1:n
            index(i,j) = index_(i,j,n);
        end
    end
end

 V = 0;
DV = zeros([dim n]);
r = zeros([dim N]);
f = zeros([1 N]);
g = zeros([1 N]);
isin = false([1 N]);
for i = 1:n-1
    k_ini = index(i,i+1);
    k_end = k_ini+n-i-1;
    r(:,k_ini:k_end) = mod(q(:,i+1:n)-q(:,i)*ones(1,n-i)+L2,L)-L2;
    for j = i+1:n
        ij = index(i,j);
%         d2 = 0;                   %
%         for k = 1:dim             % d2 = norm(r(:,ij))^2;
%             d2 = d2 + r(k,ij)^2;  % for loop faster (<) than norm
%         end                       % for < norm < '* << dot
        d2 = r(1,ij)^2 + r(2,ij)^2 + r(3,ij)^2;
        isin(ij) = (d2 < dcut2 );
        if isin(ij)
            d4 = d2^2;
            sd6 = sigma6/(d4*d2);
            V = V + sd6*(sd6-1);
            f(ij) = sd6*(1-2*sd6)/d2;
            fr = f(ij)*r(:,ij);
            DV(:,i) = DV(:,i) - fr;
            DV(:,j) = DV(:,j) + fr;
            g(ij) = sd6*(28*sd6-8)/d4;
       end
    end
end

HVv = zeros([dim n]);
DE  = zeros([dim 1]);
for i = 1:n-1
    v = DV(:,i+1:n)-DV(:,i)*ones(1,n-i);
    for j = i+1:n
        ij = index(i,j);
        if isin(ij)
            ji = j-i;
%             % DE = g(ij)*(r(:,ij)'*v(:,ji))*r(:,ij)+f(ij)*v(:,ji);
%             grv = 0;
%             for k = 1:dim
%                 grv = grv + r(k,ij)*v(k,ji);
%             end
            grv = r(1,ij)*v(1,ji) + r(2,ij)*v(2,ji) + r(3,ij)*v(3,ji);
            grv = g(ij)*grv;
%             for k = 1:dim
%                 DE(k) = grv*r(k,ij)+f(ij)*v(k,ji);
%             end
            DE(1) = grv*r(1,ij)+f(ij)*v(1,ji);
            DE(2) = grv*r(2,ij)+f(ij)*v(2,ji);
            DE(3) = grv*r(3,ij)+f(ij)*v(3,ji);
            HVv(:,i) = HVv(:,i) - DE;
            HVv(:,j) = HVv(:,j) + DE;
       end
    end
end
 V  =    4*epsilon*V;
DV  =   24*epsilon*DV;
HVv = ((24*epsilon)^2/mass)*HVv;

if iscol
    DV = DV(:);
    HVv = HVv(:);
end
end

switch order
    case 0
        % do nothing
    case 1
        V = DV;
        if nout == 2
            DV = HVv;
        end
    case 2
        V = HVv;
end
end

function k = index_(i,j,n)
k = (i-1)*n-i*(i-1)/2+(j-i);
end

function [V, DV, HV] = micLJ_(q) %#ok<DEFNU>

 V = 0;
DV = zeros([dim,n]);
HV = zeros([dim,n,dim,n]);
 I = eye(dim);
for i = 1:n-1
    r = rem(q(:,i+1:n)-q(:,i)*ones(1,n-i)+L2,L)-L2;
    for j = i+1:n
        ji = j - i;
        d = norm(r(:,ji));
        if d < dcut
            d2 = d^2;
            d4 = d2^2;
            sd6 = sigma6/(d4*d2);
            V = V + sd6*(sd6-1);
            Dv = sd6*(1-2*sd6)/d2;
            Dvr = Dv*r(:,ji);
            DV(:,i) = DV(:,i) - Dvr;
            DV(:,j) = DV(:,j) + Dvr;
            Hv = sd6*(28*sd6-8)/d4;
            Hvr = Hv*kron(r(:,ji),r(:,ji)') + Dv*I;
            Hvr = reshape(Hvr,[dim,1,dim,1]);
            HV(:,i,:,i) = HV(:,i,:,i) + Hvr;
            HV(:,i,:,j) = HV(:,i,:,j) - Hvr;
            HV(:,j,:,i) = HV(:,j,:,i) - Hvr;
            HV(:,j,:,j) = HV(:,j,:,j) + Hvr;
        end
    end
end
 V =  4*epsilon*V;
DV = 24*epsilon*DV;
HV = 24*epsilon*HV;
end
