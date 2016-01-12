function [q,p] = knd3ff(iM,fun,q0,p0,h,steps,d,k,kick1st)
%KND3FF Splitting integration method.
%   KND3FF integrates Hamilton's equations by a splitting method of kick-
%   -and-drift type (non-Hessian). It integrates iteratively and
%   alternately each equation. Properties: 4th order (effective), 3.5-
%   -stage, explicit, self-adjoint, non-symplectic, fast processing. For
%   more details, see [1-2].
%
%   Call:
%
%      [q,p] = KND3FF(iM,fun,q0,p0,h,steps,d,k,kick1st)
%
%   Inputs:
%
%      iM, fun, q0, p0, h, steps, kick1st as in ROWLANDS.
%
%      d   Main drift coefficient, a real number. The others are computed
%          upon execution to achieve effectiveness.
%      k   Main kick coefficient, otherwise as input 'd' above.
%
%   Outputs:
%
%      q, p as in ROWLANDS.
%
%   Example:
%
%      >> [q,p] = knd3ff(1,@(q)q,2,2,0.1,35,1/2,1/2);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with a drift1st integrator and plots the resulting trajectory in
%      phase space.
%
%   See also KND, KND3MS, HKND2FF, ROWLANDS.
%
%   [1] M.A. Lopez-Marcos, J.M. Sanz-Serna, R.D. Skeel: "Cheap enhancement
%   of symplectic integrators", Numerical Analysis (1995, Dundee), Pitman
%   Res. Notes Math. 344 (1996)
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
%   $Revision: 0.1.0.01 $  $Date: 2016/01/11 17:34:52 $

% parameters
if exist('kick1st','var') && kick1st
     kick = [k 1/2-k 1/2-k k];
    drift = [1/2-d 2*d 1/2-d];
%      l   = 1/24-(1/2-k)^2*(1/2-d);
      l   = (1/2-d)*(1/2-k)*(k-d-1)/2+1/16;
else % drift1st
%     drift = [1/2-d d d 1/2-d];
%      kick = [  k  1-2*k  k  ];
%       l   = d*k*(k-d-1)/2+1/16;
%       l   = -l;
    drift = [d 1/2-d 1/2-d d];
     kick = [1/2-k 2*k 1/2-k];
%      l   = (1/2-d)^2*(1/2-k)-1/24;
      l   = (1/2-d)*(1/2-k)*(k-d+1)/2-1/16;
end
hm = h/2;

% pre-processing
fp = fun(q0+hm*iM*p0);
fm = fun(q0-hm*iM*p0);

Q0 = q0 + h^2*l*iM*(fp+fm)/2;
P0 = p0 -   h*l*(fp-fm);

% integration
[q,p] = knd(iM,fun,Q0,P0,h,steps+1,drift,kick);

% post-processing
q = q(:,1:steps)+l*(q(:,2:steps+1)-2*q(:,1:steps)+[Q0,q(:,1:steps-1)]);
p = p(:,1:steps)-l*(p(:,2:steps+1)-2*p(:,1:steps)+[P0,p(:,1:steps-1)]);
end
