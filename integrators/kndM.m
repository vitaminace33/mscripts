function M = kndM(h,drift,kick,hessian,kick1st)
%KNDM Amplification matrix for splitting methods.
%   Consider the basic Hamiltonian model system, the harmonic oscillator,
%
%      H(q,p) = 1/2*p^2 + 1/2*q^2 .
%
%   KNDM computes the amplification matrix of this system for r-stages
%   splitting methods of (Hessian)-kick-and-drift type, useful for the
%   method's analysis. For more details, see [1-2].
%
%   Call:
%
%      [q,p] = KNDM(h,drift,kick,hessian,kick1st)
%
%   Inputs:
%
%      h, drift, kick, hessian, kick1st as in KND.
%
%   Outputs:
%
%      M   The amplification matrix, a 2x2 matrix.
%
%   Example:
%
%      >> syms h positive; syms c real;
%      >> M = kndM(h,[1/2 1/2],1,c);
%
%      The previous code computes the amplification Matrix for a position-
%      -Rowlands-like integrator, which coincides with Rolands' method
%      itself for c=-1/24, with Verlet's integrator for c=0, with a method
%      of 2nd-order/4th-order in phase/energy for c=-3/8.
%
%   See also KND.
%
%   [1] M.A. Lopez-Marcos, J.M. Sanz-Serna, R.D. Skeel: "An explicit
%   symplectic integrator with maximal stability interval", Numerical
%   Analysis (A.R. Mitchell 75th B-day), World Sci. Publ. (1996)
%
%   [2] J.M. Sanz-Serna, M.P. Calvo: "Numerical Hamiltonian problems",
%   Chapman & Hall (1994)
%
%   URLs: https://github.com/vitaminace33
%         https://google.com/+cedricmartinezcampos
%
%   Copyleft 2014-2016 Cedric M. Campos, Dr.
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
%   $Revision: 0.2.0.01 $  $Date: 2016/01/11 17:34:56 $
hd = h*drift(:);
if nargin>3
    kick = kick+2*h^2*hessian;
end
hk = h*kick(:);
m = size(hd,1);
n = size(hk,1);
switch n-m
    case 1
        kick1st  = true;
        kicklast = true;
        l = m;
    case -1
        kick1st  = false;
        kicklast = false;
        l = n;
    case 0
        kicklast = ~kick1st;
        l = n;
    otherwise
        error('integrators:kndM:inconsistency', ...
              ['Inconsistent kick-drift parameters', ...
               '; their lengths must differ at most in 1.'])
end
M = eye(2);
if kick1st % momemtum kick 1st
    for i = 1:l
        M = kickM(hk(i))*M;
        M = driftM(hd(i))*M;
    end
    if kicklast % momemtum kick last
        M = kickM(hk(l+1))*M;
    end
else % position drift 1st
    for i = 1:l
        M = driftM(hd(i))*M;
        M = kickM(hk(i))*M;
    end
    if ~kicklast % position drift last
        M = driftM(hd(l+1))*M;
    end
end
end

function D = driftM(hd)
D = [1 hd;0 1];
end

function K = kickM(hk)
K = [1 0;-hk 1];
end
