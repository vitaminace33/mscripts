function [q,p] = verlet(iM,fun,q0,p0,h,steps,f0)
%VERLET Stormer-Verlet/leapfrog splitting integration method.
%   VERLET integrates Hamilton's equations by a splitting method of kick-
%   -and-drift type (non-Hessian). It integrates iteratively and
%   alternately each equation. Properties: 2nd order, 1.5-stage, explicit,
%   self-adjoint, symplectic. For more details, see [1].
%
%   Call:
%
%      [q,p] = VERLET(iM,fun,q0,p0,h,steps,f0)
%
%   Inputs:
%
%      iM, fun, q0, p0, h, steps, f0 as in KND.
%
%   Outputs:
%
%      q, p as in KND.
%
%   Example:
%
%      >> [q,p] = verlet(1,@(q)q,2,2,0.1,35);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with the position-Verlet integrator and plots the resulting
%      trajectory in phase space.
%
%   See also KND.
%
%   [1] J.M. Sanz-Serna, M.P. Calvo: "Numerical Hamiltonian problems",
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
%   $Revision: 0.1.0.01 $  $Date: 2016/01/11 17:34:48 $
hm = h/2;
switch length(steps)
    case 1
        stepsps = 1;
        samples = steps;
    case 2
        stepsps = steps(1);
        samples = steps(2);
    otherwise
        error('integrators:verlet:wrongsteps', ...
              ['Wrong steps/samples parameters', ...
               '; tsteps must have length 1 or 2.'])
end
if exist('f0','var')
    kick1st  = true;
else
    kick1st = false;
end
if samples == 1
    q = q0; p = p0;
    if kick1st
        if isempty(f0)
            f = fun(q0);
        else
            f = f0;
        end
        p = p - hm*f;
        for i = 1:steps-1
            q = q + h*iM*p;
            p = p - h*fun(q);
        end
        q = q + h*iM*p;
        p = p - hm*fun(q);
    else
        q = q + hm*iM*p;
        for i = 1:steps-1
            p = p - h*fun(q);
            q = q + h*iM*p;
        end
        p = p - h*fun(q);
        q = q + hm*iM*p;
    end
else
    d = size(q0,1);
    q = zeros([d samples]);
    p = zeros([d samples]);
    q1 = q0; p1 = p0;
    if kick1st
        if isempty(f0)
            f1 = fun(q0);
        else
            f1 = f0;
        end
        for i = 1:samples
            p1 = p1 - hm*f1;
            for j = 1:stepsps-1
                q1 = q1 + h*iM*p1;
                p1 = p1 - h*fun(q1);
            end
            q1 = q1 + h*iM*p1;
            f1 = fun(q1);
            p1 = p1 - hm*f1;
            q(:,i) = q1; p(:,i) = p1;
        end
    else
        for i = 1:samples
            q1 = q1 + hm*iM*p1;
            for j = 1:stepsps-1
                p1 = p1 - h*fun(q1);
                q1 = q1 + h*iM*p1;
            end
            p1 = p1 - h*fun(q1);
            q1 = q1 + hm*iM*p1;
            q(:,i) = q1; p(:,i) = p1;
        end
    end
end
end
