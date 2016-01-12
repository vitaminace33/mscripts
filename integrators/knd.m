function [q,p] = knd(iM,fun,q0,p0,h,steps,drift,kick,hessian,f0,g0)
%KND Splitting method to integrate mechanical systems.
%   Consider an n-dimensional mechanical system whose total energy is given
%   by the Hamiltonian function
%
%      H(q,p) = T(p) + V(q) ,   with T(q,p) = 1/2*p'*iM*p ,
%
%   where (q,p) are respectively the current configuration and momentum,
%   T(p) is the kinetic energy (iM being the inverse of the mass matrix)
%   and V(q) the potential energy. The natural motions of the system are
%   characterized by Hamilton's equations
%
%      dq/dt = D_pH = iM*p ,   dp/dt = -D_qH = -DV .
%
%   KND integrates these by a splitting method of (Hessian)-kick-and-drift
%   type. It integrates iteratively and alternately each equation using,
%   optionally, the Hessian of the potential.
%
%   The order of the method, as well as other properties (simplecticity,
%   energy preservation), depend on the integration coefficients. For more
%   details, see [1].
%
%   Call:
%
%      [q,p] = KND(iM,fun,q0,p0,h,steps,drift,kick,hessian,f0,g0)
%
%   Inputs:
%
%         iM   Inverse mass matrix. Square matrix of dimension n. If
%              scalar, set to the only eigenvalue to lower load.
%        fun   Potential gradient. Function handler that applied to a
%              configuration q, outputs DV(q), a column-vector of length n.
%              If a Hessian method is chosen, it must output as second
%              element HV*iM*DV, also a column-vector of length n.
%         q0   Initial configuration. Column-vector of length n.
%         p0   Initial momentum. Column-vector of length n.
%          h   Time step. Positive scalar.
%      steps   Number of steps. It can be either an integer with the total
%              number of steps or the couple of integers [stepsps samples].
%                 stepsps   Steps per sample. Default value: 1.
%                 samples   Number of outputted samples (see below).
%                           Default value: "steps".
%              The total number of steps is therefore stepsps*samples.
%      drift   Drift coefficients. Vector of unitary sum.
%       kick   Kick coefficients. Vector of unitary sum. Its length differs
%              (at most) in 1 with respect to the length of "drift".
%    hessian   Hessian coefficients (optional). Vector of same length than
%              "kick". Set it to the empty array [] to force a non-Hessian
%              method.
%         f0   Gradient at q0 (optional), see "fun" above. If set, it
%              implies that "hessian" is set as-well. It is only used in
%              kick-first methods.
%         g0   Hessian-gradient at q0 (optional), see "fun" above. If set,
%              it implies that "hessian" and "f0" are set as-well. It is
%              only used in Hessian-kick-first methods.
%
%   Outputs:
%
%          q   Computed configuration samples. Array of size [n samples].
%          p   Computed momentum samples. Array of size [n samples].
%
%   Note: let n=length(kick), m=length(drift), r=(n+m)/2, then the method
%   is an r-stages splitting one such that
%      * if n=m+1, it is of kick-1st type (as well of kick-last);
%      * if m=n+1, it is of drift-1st type (as well of drift-last);
%      * if n=m, the kick-1st type (and drift-last) is implied setting up
%        the variable 'f0', otherwise the drift-1st type (kick-last) is
%        considered.
%
%   Example:
%
%      >> [q,p] = knd(1,@(q)q,2,2,0.1,35,[1/2 1/2],1);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with the position-Verlet integrator and plots the resulting
%      trajectory in phase space.
%
%   See also VERLET, ROWLANDS.
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
%   $Revision: 0.2.0.01 $  $Date: 2016/01/11 17:34:47 $
switch length(steps)
    case 1
        stepsps = 1;
        samples = steps;
    case 2
        stepsps = steps(1);
        samples = steps(2);
    otherwise
        error('integrators:knd:wrongsteps', ...
              ['Wrong steps/samples parameters', ...
               '; ''steps'' must have length 1 or 2.'])
end
hd = h*drift(:);
hk = h*kick(:);
m = size(hd,1);
n = size(hk,1);
switch n-m
    case 1
        kick1st  = true;
        kicklast = true;
        l = m;
%         hk1n = h*(kick(1)+kick(n));
        if ~exist('f0','var')
            f0 = [];
        end
    case -1
        kick1st  = false;
        kicklast = false;
        l = n;
%         hd1m = h*(drift(1)+drift(m));
    case 0
        if exist('f0','var')
            kick1st  = true;
            kicklast = false;
        else
            kick1st  = false;
            kicklast = true;
        end
        l = m; % = n
    otherwise
        error('integrators:knd:inconsistency', ...
              ['Inconsistent kick-drift parameters', ...
               '; their lengths must differ at most in 1.'])
end

if exist('hessian','var') && ~isempty(hessian)
    ishessian = true;
    h3c2 = 2*h^3*hessian(:);
    if ~exist('g0','var')
        g0 = [];
    end
    if size(h3c2,1) ~= n
        error('integrators:hknd:inconsistency', ...
              ['Inconsistent hessian-kick parameters', ...
               '; their lengths must coincide.'])
    end
else
    ishessian = false;
end
dim = size(q0,1);
q = zeros([dim samples]);
p = zeros([dim samples]);
q1 = q0;
p1 = p0;
if ishessian % hessian method
    if kick1st % momemtum kick 1st
        if isempty(f0) || isempty(g0)
            [f1,g1] = fun(q0);
        else
            f1 = f0;
            g1 = g0;
        end
        if kicklast % momemtum kick last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        p1 = p1 - (hk(k)*f1 + h3c2(k)*g1);
                        q1 = q1 + hd(k)*iM*p1;
                        [f1,g1] = fun(q1);
                    end
                    p1 = p1 - (hk(n)*f1 + h3c2(n)*g1);
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        else % position drift last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        p1 = p1 - (hk(k)*f1 + h3c2(k)*g1);
                        q1 = q1 + hd(k)*iM*p1;
                        [f1,g1] = fun(q1);
                    end
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        end
    else % position drift 1st
        if kicklast % momemtum kick last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        q1 = q1 + hd(k)*iM*p1;
                        [f1,g1] = fun(q1);
                        p1 = p1 - (hk(k)*f1 + h3c2(k)*g1);
                    end
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        else % position drift last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        q1 = q1 + hd(k)*iM*p1;
                        [f1,g1] = fun(q1);
                        p1 = p1 - (hk(k)*f1 + h3c2(k)*g1);
                    end
                    q1 = q1 + hd(m)*iM*p1;
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        end
    end
else % regular method (non-hessian)
    if kick1st % momemtum kick 1st
        if isempty(f0)
            f1 = fun(q0);
        else
            f1 = f0;
        end
        if kicklast % momemtum kick last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        p1 = p1 - hk(k)*f1;
                        q1 = q1 + hd(k)*iM*p1;
                        f1 = fun(q1);
                    end
                    p1 = p1 - hk(n)*f1;
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        else % position drift last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        p1 = p1 - hk(k)*f1;
                        q1 = q1 + hd(k)*iM*p1;
                        f1 = fun(q1);
                    end
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        end
    else % position drift 1st
        if kicklast % momemtum kick last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        q1 = q1 + hd(k)*iM*p1;
                        f1 = fun(q1);
                        p1 = p1 - hk(k)*f1;
                    end
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        else % position drift last
            for i = 1:samples
                for j = 1:stepsps
                    for k = 1:l
                        q1 = q1 + hd(k)*iM*p1;
                        f1 = fun(q1);
                        p1 = p1 - hk(k)*f1;
                    end
                    q1 = q1 + hd(m)*iM*p1;
                end
                q(:,i) = q1;
                p(:,i) = p1;
            end
        end
    end
end
end
