function [q, p, accepted, N] = xhmc(funs, q0, options)
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
%   $Revision: 0.1.0.00 $  $Date: 2016/01/11 18:01:02 $

% default options
%      N = 1000; % Markov chain length
%   burn = 100;  % MCMC burn-in length
%   beta = 1;    % Boltzmann inverse temperature
%    psi = pi/3; % Horowitz's angle, psi = pi/2 gives HMC
%  extra = 3;    % Extra chances, extra = 0 gives (G)HMC
% integrator = @verlet;
%      h = .01;  % Time step integration
%  steps = 10;   % Integration steps
%  shift = .1;   % ±shift% for t.s.i. above (see also `h' below)
% MaxInt = inf;  % Maximum number of integrations

eo = isempty(options);
if  eo || ~isfield(options,'N') || isempty(options.N)
    options.N = 1000;
end
if  eo || ~isfield(options,'burn') || isempty(options.burn)
    options.burn = 100;
end
if  eo || ~isfield(options,'beta') || isempty(options.beta)
    options.beta = 1;
end
if  eo || ~isfield(options,'psi') || isempty(options.psi)
    options.psi = pi/4;
end
if  eo || ~isfield(options,'extra') || isempty(options.extra)
    options.extra = 3;
end
if  eo || ~isfield(options,'integrator') || isempty(options.integrator)
    options.integrator = @verlet;
end
if  eo || ~isfield(options,'h') || isempty(options.h)
    options.h = .01;
end
if  eo || ~isfield(options,'steps') || isempty(options.steps)
    options.steps = 10;
end
if  eo || ~isfield(options,'shift') || isempty(options.shift)
    options.shift = .1;
end
if  eo || ~isfield(options,'MaxInt') || isempty(options.MaxInt)
    options.MaxInt = inf;
end
if  eo || ~isfield(options,'iM') || isempty(options.iM)
    options.iM = 1;
end

integrator = options.integrator;
     iM = options.iM;
 MaxInt = options.MaxInt;
   beta = options.beta;
    psi = options.psi;
      N = floor(abs(options.N));
   burn = floor(abs(options.burn));
      h = options.h;
  steps = options.steps;
  shift = options.shift;
  tries = options.extra+1;

cospsi = cos(psi);
sinpsi = sin(psi);

% dimension
dim = size(q0);
if length(dim) > 2 || dim(2)~=1
    error('Q0 must be a column vector.')
end
dim = dim(1);

% potential energy (and gradient)
if isa(funs,'function_handle')
     V = @(q)funs(q,0);
    DV = @(q)funs(q,1);
elseif isa(funs,'cell') && length(funs)==2
     V = funs{1};
    DV = funs{2};
else
    error('FUNS must be a function handle or a cell of (2) function handles.')
end

% initial burn-in phase (recursive)
if burn > 0
    options.burn = 0;
    options.N = burn;
    options.MaxInt = Inf;
    q = xhmc(funs,q0,options);
    q0 = q(:,burn);
end

irbeta = 1/sqrt(beta);

q = zeros([dim N]);
p = zeros([dim N]);

q1 = q0;
p1 = zeros([dim 1]);
V1 = V(q0);
accepted = zeros(1,tries);
i = 0;
NumInt = 0;
while i < N && NumInt < MaxInt
    i = i+1;
    q0 = q1;
    p0 = p1;
    V0 = V1;
    p0 = cospsi*p0 + sinpsi*irbeta*randn([dim 1]);
    E0 = 1/2*(p0'*iM*p0) + V0;
    u = rand;
    trynum = 0;
    acceptance = 0;
    hplus = h*(1-shift+2*shift*rand);
    q1 = q0;
    p1 = p0;
    while acceptance < u && trynum < tries
        trynum = trynum + 1;
        [q1, p1] = integrator(iM,DV,q1,p1,hplus,steps);
        V1 = V(q1);
        E1 = 1/2*(p1'*iM*p1) + V1;
        acceptance = min(1,max(acceptance,exp(-beta*(E1-E0))));
        E0 = E1;
    end
    if acceptance >= u
        accepted(trynum) = accepted(trynum) +1;
    else
        q1 =  q0;
        p1 = -p0;
        V1 =  V0;
    end
    q(:,i) = q1;
    p(:,i) = p1;
    NumInt = NumInt + trynum;
end
N = i;
end

function [q1,p1] = verlet(iM,DV,q0,p0,h,steps)
q1 = q0;
pm = p0 - h/2*DV(q0); % DV(q0) could reuse DV(q1) from past integration
for i = 1:steps-1
    q1 = q1 + h*iM*pm;
    pm = pm - h*DV(q1);
end
q1 = q1 + h*iM*pm;
p1 = pm - h/2*DV(q1); % DV(q1) could be salvaged at next integration
end
