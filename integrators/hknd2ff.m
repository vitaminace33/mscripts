function [q,p] = hknd2ff(iM,fun,q0,p0,h,steps,params)
%HKND2FF Splitting integration method.
%   HKND2FF integrates Hamilton's equations by a splitting method of
%   Hessian-kick-and-drift type. It integrates iteratively and alternately
%   each equation with the aid of the Hessian. Properties: 4th order
%   (effective), 2.5-stage, explicit, self-adjoint, non-symplectic, fast
%   processing. For more details, see [1-3].
%
%   Call:
%
%      [q,p] = HKND2FF(iM,fun,q0,p0,h,steps,params)
%
%   Inputs:
%
%      iM, fun, q0, p0, h, steps as in KND/HKND2SA.
%
%       steps   Number of total steps and samples, i.e. 1 sample per step.
%               (Differs from KND/HKND2SA due to processing)
%      params   Vector of Hessian-kick-and-drift parameters. It can either
%               be [d] or [k c].
%
%               d   Main drift coefficient, a real number. The other
%                   coefficients are computed upon execution to achieve
%                   self-adjointness and effectiveness.
%               k   Main kick coefficient, otherwise as input 'd' above.
%               c   Main Hessian coefficient, otherwise as input 'd' above.
%
%   Outputs:
%
%      q, p as in KND/HKND2SA.
%
%      >> [q,p] = hknd2sa(1,@harmosc,2,2,0.1,35,1/4);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with a double-concatenation of the position-Rowlands integrator
%      (with processing) and plots the resulting trajectory in phase
%      space.
%
%   See also KND, KND3FF, HKND2SA, ROWLANDS.
%
%   [1] J.M. Sanz-Serna, M.P. Calvo: "Numerical Hamiltonian problems",
%   Chapman & Hall (1994)
%
%   [2] M.A. Lopez-Marcos, J.M. Sanz-Serna, R.D. Skeel: "Explicit
%   Symplectic Integrators Using Hessian--Vector Products", SIAM J. Sci.
%   Comput. 18 (1997)
%
%   [3] J.M. Sanz-Serna, M.P. Calvo: "Numerical Hamiltonian problems",
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
%   $Revision: 0.1.0.01 $  $Date: 2016/01/11 17:34:58 $
switch length(params)
    case 1
        kick1st = false;
        % a = params;
        a = params;
        c = -( (a-1/4)^2+1/48 )/4; % check the leading sign!!!
        l = ( (a-1/2)^2-1/12 )/2;
    case 2
        kick1st = true;
        % [b,co] = deal(params);
        b  = params(1);
        co = params(2);
        ci = -( 2*co + (b-1/4)^2/2 + 1/96 ); % check the leading sign!!!
        l  = (b-1/6)/4;
    otherwise
        error('integrators:hknd2ff:inconsistency', ...
              ['Inconsistent Hessian-kick-drift parameter', ...
               '; ''parameter'' must have length 1 or 2.'])
end
hm = h/2;

% pre-processing
fp = fun(q0+hm*iM*p0);
fm = fun(q0-hm*iM*p0);

Q0 = q0 + h^2*l*iM*(fp+fm)/2;
P0 = p0 -   h*l*iM*(fp-fm);

% integration
if kick1st
  [q,p] = hknd2sa(iM,fun,Q0,P0,h,steps+1,[b,co,ci]);
else
  [q,p] = hknd2sa(iM,fun,Q0,P0,h,steps+1,[a,c]);
end

% post-processing
q = q(:,1:steps)+l*(q(:,2:steps+1)-2*q(:,1:steps)+[Q0,q(:,1:steps-1)]);
p = p(:,1:steps)-l*(p(:,2:steps+1)-2*p(:,1:steps)+[P0,p(:,1:steps-1)]);
end
