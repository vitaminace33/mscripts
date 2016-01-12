function [q,p] = knd3nrg(iM,fun,q0,p0,h,steps)
%KND3NRG Splitting integration method.
%   KND3NRG integrates Hamilton's equations by a splitting method of kick-
%   -and-drift type (non-Hessian). It integrates iteratively and
%   alternately each equation. Properties: 2nd order (phase),  4th order
%   (energy), 3.5-stage, explicit, self-adjoint, symplectic. For more
%   details, see [1-2].
%
%   Call:
%
%      [q,p] = KND3NRG(iM,fun,q0,p0,h,steps)
%
%   Inputs:
%
%      iM, fun, q0, p0, h, steps (as in KND)
%
%   Outputs:
%
%      q, p as in KND.
%
%   Example:
%
%      >> [q,p] = knd3nrg(1,@(q)q,2,2,0.1,35);
%      >> plot(q,p)
%
%      The previous code integrates a simple harmonic oscillator system
%      with a drift1st integrator and plots the resulting trajectory in
%      phase space.
%
%   See also KND, KND3NRG.
%
%   [1] C. Predescu, et al.: "Computationally efficient molecular dynamics
%   integrators with improved sampling accuracy", Mol. Phys. 110 (2012)
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
%   $Revision: 0.1.0.01 $  $Date: 2016/01/11 17:34:54 $

% parameters
% 4th order (energy)
a = 0.108991425403425322;
b = 0.290485609075128726;
% % 6th order (energy)
% a = 0.099031132097580873;
% b = 0.277479066043685596;
% % drift1st
% drift = [a 1/2-a 1/2-a a];
%  kick = [ b   1-2*b   b ];
% kick1st
 kick = [a 1/2-a 1/2-a a];
drift = [ b   1-2*b   b ];

% integration
[q,p] = knd(iM,fun,q0,p0,h,steps,drift,kick);
end
