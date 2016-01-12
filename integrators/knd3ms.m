function [q,p] = knd3ms(iM,fun,q0,p0,h,steps,kick1st)
%KND3MS Splitting integration method.
%   KND3MS integrates Hamilton's equations by a splitting method of kick-
%   -and-drift type (non-Hessian). It integrates iteratively and
%   alternately each equation. Properties: 4th order (effective), 3.5-
%   -stage, self-adjoint, non-symplectic, fast processing, maximal
%   stability. For more details, see [1-4].
%
%   Call:
%
%      [q,p] = KND3MS(iM,fun,q0,p0,h,steps,kick1st)
%
%   Inputs:
%
%      iM, fun, q0, p0, h, steps, kick1st as in ROWLANDS.
%
%   Outputs:
%
%      q, p as in ROWLANDS.
%
%   Example:
%
%      >> [q,p] = knd3ms(1,@(q)q,2,2,0.1,35);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with a drift1st integrator and plots the resulting trajectory in
%      phase space.
%
%   See also KND, KND3FF, HKND2FF, ROWLANDS.
%
%   [1] M.A. Lopez-Marcos, J.M. Sanz-Serna, R.D. Skeel: "Cheap enhancement
%   of symplectic integrators", Numerical Analysis (1995, Dundee), Pitman
%   Res. Notes Math. 344 (1996)
%
%   [2] M.A. Lopez-Marcos, J.M. Sanz-Serna, R.D. Skeel: "An explicit
%   symplectic integrator with maximal stability interval", Numerical
%   Analysis (A.R. Mitchell 75th B-day), World Sci. Publ. (1996)
%
%   [3] J.M. Sanz-Serna, M.P. Calvo: "Numerical Hamiltonian problems",
%   Chapman & Hall (1994)
%
%   [4] H. Yoshida: "Construction of higher order symplectic integrators",
%   Phys. Lett. A 150 (1990)
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
%   $Revision: 0.1.0.01 $  $Date: 2016/01/11 17:34:50 $

% parameters
cr2 = 2^(1/3);
  a = (2+cr2+1/cr2)/6;
  b = 1/2-a;
if exist('kick1st','var') && kick1st
    % kick 1st
     kick = [a b b a];
    drift = [b 2*a b];
      l = (1+cr2)/48;
else
    % drift 1st
    drift = [a b b a];
     kick = [b 2*a b];
      l = -(1+cr2)/48;
end
 hm = h/2;

% pre-processing
fp = fun(q0+hm*iM*p0);
fm = fun(q0-hm*iM*p0);

Q0 = q0 + h^2*l*iM*(fp+fm)/2;
P0 = p0 -   h*l*iM*(fp-fm);

% integration
[q,p] = knd(iM,fun,Q0,P0,h,steps+1,drift,kick);

% post-processing
q = q(:,1:steps)+l*(q(:,2:steps+1)-2*q(:,1:steps)+[Q0,q(:,1:steps-1)]);
p = p(:,1:steps)-l*(p(:,2:steps+1)-2*p(:,1:steps)+[P0,p(:,1:steps-1)]);
end
