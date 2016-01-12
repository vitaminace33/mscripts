function [q,p] = hknd2sa(iM,fun,q0,p0,h,steps,params,f0,g0)
%HKND2SA Splitting integration method.
%   HKND2SA integrates Hamilton's equations by a splitting method of
%   Hessian-kick-and-drift type. It integrates iteratively and alternately
%   each equation with the aid of the Hessian. Properties: 2nd order
%   2.5-stage, explicit, self-adjoint, symplectic. For more details, see
%   [1-3].
%
%   Call:
%
%      [q,p] = HKND2SA(iM,fun,q0,p0,h,steps,params,f0,g0)
%
%   Inputs:
%
%      iM, fun, q0, p0, h, steps, f0, g0 as in KND.
%
%      params   Vector of Hessian-kick-and-drift parameters. It can either
%               be [d c] or [k co ci].
%
%                d   Main drift coefficient, a real number. The others are
%                    computed upon execution to achieve self-adjointness.
%                k   Main kick coefficient, otherwise as input 'd' above.
%               c?   Hessian coefficients, real numbers. Respectively, the
%                    main and only one, the main and outer one, the
%                    secondary and inner one.
%
%   Outputs:
%
%      q, p as in KND.
%
%   Note: the Hessian-vector product is used iff the corresponding Hessian
%   component is non-null.
%
%   Example:
%
%      >> [q,p] = hknd2sa(1,@harmosc,2,2,0.1,35,[1/4 -1/192]);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with a double-concatenation of the position-Rowlands integrator
%      (without processing) and plots the resulting trajectory in phase
%      space.
%
%   See also KND, KND3FF, HKND2FF, ROWLANDS.
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
%   $Revision: 0.1.0.01 $  $Date: 2016/01/11 17:34:57 $
switch length(steps)
    case 1
        stepsps = 1;
        samples = steps;
    case 2
        stepsps = steps(1);
        samples = steps(2);
    otherwise
        error('integrators:hknd2sa:wrongsteps', ...
              ['Wrong steps/samples parameters', ...
               '; ''steps'' must have length 1 or 2.'])
end
switch length(params)
    case 2
        kick1st = false;
        % [d,c] = deal(params);
        hm  = h/2;
        hdo =       params(1)*h;
        hdi = (1-2*params(1))*h;
        h3c =     2*params(2)*h^3;
    case 3
        kick1st = true;
        % [k,co,ci] = deal(params);
        hm   = h/2;
        hko  =       params(1)*h;
        hki  = (1-2*params(1))*h;
        h3co =     2*params(2)*h^3;
        h3ci =     2*params(3)*h^3;
        if ~exist('f0','var')
            f0 = [];
        end
        if ~exist('g0','var')
            g0 = [];
        end
    otherwise
        error('integrators:hknd2sa:inconsistency', ...
              ['Inconsistent Hessian-kick-drift parameter', ...
               '; ''parameter'' must have length 2 or 3.'])
end

dim = size(q0,1);
q = zeros([dim samples]);
p = zeros([dim samples]);
q1 = q0;
p1 = p0;

if kick1st
    if params(3) == 0
        if isempty(f0) || isempty(g0)
            [f1,g1] = fun(q0);
        else
            f1 = f0;
            g1 = g0;
        end
        for i = 1:samples
            for j = 1:stepsps
                p1 = p1 - (hko*f1 + h3co*g1);
                q1 = q1 + hm*iM*p1;
                f1 = fun(q1);
                p1 = p1 - hki*f1;
                q1 = q1 + hm*iM*p1;
                [f1,g1] = fun(q1);
                p1 = p1 - (hko*f1 + h3co*g1);
            end
            q(:,i) = q1;
            p(:,i) = p1;
        end
    elseif params(2) == 0
        if isempty(f0)
            f1 = fun(q0);
        else
            f1 = f0;
        end
        for i = 1:samples
            for j = 1:stepsps
                p1 = p1 - hko*f1;
                q1 = q1 + hm*iM*p1;
                [f1,g1] = fun(q1);
                p1 = p1 - (hki*f1 + h3ci*g1);
                q1 = q1 + hm*iM*p1;
                f1 = fun(q1);
                p1 = p1 - hko*f1;
            end
            q(:,i) = q1;
            p(:,i) = p1;
        end
    else
        if isempty(f0) || isempty(g0)
            [f1,g1] = fun(q0);
        else
            f1 = f0;
            g1 = g0;
        end
        for i = 1:samples
            for j = 1:stepsps
                p1 = p1 - (hko*f1 + h3co*g1);
                q1 = q1 + hm*iM*p1;
                [f1,g1] = fun(q1);
                p1 = p1 - (hki*f1 + h3ci*g1);
                q1 = q1 + hm*iM*p1;
                [f1,g1] = fun(q1);
                p1 = p1 - (hko*f1 + h3co*g1);
            end
            q(:,i) = q1;
            p(:,i) = p1;
        end
    end
else
    for i = 1:samples
        for j = 1:stepsps
            q1 = q1 + hdo*iM*p1;
            [f1,g1] = fun(q1);
            p1 = p1 - (hm*f1 + h3c*g1);
            q1 = q1 + hdi*iM*p1;
            [f1,g1] = fun(q1);
            p1 = p1 - (hm*f1 + h3c*g1);
            q1 = q1 + hdo*iM*p1;
        end
        q(:,i) = q1;
        p(:,i) = p1;
    end
end
end
