function [q,p] = rowlands(iM,fun,q0,p0,h,steps,kick1st)
%ROWLANDS Splitting integration method.
%   ROWLANDS integrates Hamilton's equations by a splitting method of
%   Hessian-kick-and-drift type. It integrates iteratively and alternately
%   each equation with the aid of the Hessian. Properties: 4th order
%   (effective), 1.5-stage, explicit, self-adjoint, non-symplectic, fast
%   processing. For more details, see [1-3].
%
%   Call:
%
%      [q,p] = ROWLANDS(iM,fun,q0,p0,h,steps,kick1st)
%
%   Inputs:
%
%      iM, fun, q0, p0, h as in KND.
%
%        steps   Number of total steps and samples, i.e. 1 sample per step.
%                (Differs from KND due to processing)
%      kick1st   Boolean flag for method's type (optional). Default value:
%                false. (Differs from KND due to processing)
%
%   Outputs:
%
%      q, p as in KND.
%
%   Example:
%
%      >> [q,p] = rowlands(1,@harmosc,2,2,0.1,35);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with the position-Rowlands integrator and plots the resulting
%      trajectory in phase space. Above we assume:
%
%      function [f,g] = harmosc(q)
%          f = q;    g = q;
%      end
%
%   See also KND.
%
%   [1] M.A. Lopez-Marcos, J.M. Sanz-Serna, R.D. Skeel: "Cheap enhancement
%   of symplectic integrators", Numerical Analysis (1995, Dundee), Pitman
%   Res. Notes Math. 344 (1996)
%
%   [2] G. Rowlands: "A numerical algorithm for Hamiltonian systems", J.
%   Comput. Phys. 97 (1991)
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
%   $Revision: 0.1.0.01 $  $Date: 2016/01/11 17:34:49 $
if exist('kick1st','var') && kick1st
    l = 1/12;
else % drift1st
    kick1st = false;
    l  = -1/24;
end
hm = h/2;
h2c = 2*h^2/24;

% pre-processing
fp = fun(q0+hm*iM*p0);
fm = fun(q0-hm*iM*p0);

Q0 = q0 + h^2*l*iM*(fp+fm)/2;
P0 = p0 -   h*l*iM*(fp-fm);

% integration
% if kick1st
%     [q,p] = knd(iM,fun,Q0,P0,h,steps+1,1,[1,1]/2,-[1,1]/24/2);
% else
%     [q,p] = knd(iM,fun,Q0,P0,h,steps+1,[1,1]/2,1,-1/24);
% end
% q = [Q0,q]; p = [P0,p];
d = size(q0,1);
q = zeros([d steps+2]);
p = zeros([d steps+2]);
q(:,1) = Q0;
p(:,1) = P0;
if kick1st
    [fj,gj] = fun(q(:,1));
    i = 1;
    for j = 2:steps+2
        pm     = p(:,i) - hm*(fj-h2c*gj);
        q(:,j) = q(:,i) + h*iM*pm;
        [fj,gj] = fun(q(:,j));
        p(:,j) = pm     - hm*(fj-h2c*gj);
        i = j;
    end
else
    i = 1;
    for j = 2:steps+2
        qm     = q(:,i) + hm*iM*p(:,i);
        [fm,gm] = fun(qm);
        p(:,j) = p(:,i) - h*(fm-h2c*gm);
        q(:,j) = qm     + hm*iM*p(:,j);
        i = j;
    end
end

% post-processing
q = q(:,2:steps+1)+l*(q(:,3:steps+2)-2*q(:,2:steps+1)+q(:,1:steps));
p = p(:,2:steps+1)-l*(p(:,3:steps+2)-2*p(:,2:steps+1)+p(:,1:steps));
end
