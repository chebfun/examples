%% Linear `exp` initial-value problem
% Tom Maerz, October 2010

%%
% (Chebfun example ode-linear/LinExpIVP.m)
% [Tags: #linearODE, #IVP]

%%
% This is an elementary example to illustrate how one might use Chebfun to
% solve a very simple ODE initial-value problem. We take the scalar test
% problem
%
% $$ u' - \lambda u = 0  ,~~~    u(0) = 1,~  \lambda = -10000 $$
%
% on the interval $[0,.005]$.  The solution is $\exp(\lambda x)$.
d = [0,.005];                       % domain
x = chebfun('x', d);                % x variable
L = chebop(d);                      % operator
lambda = -10000;                    % specifying parameter lambda
L.op = @(u) diff(u,1) - lambda*u;   % linear operator defining the ODE
L.lbc = @(u) u-1;                   % imposing Dirichlet boundary condition
u = L\0;                            % solve the problem
plot(u,'linewidth',1.6)             % plot the solution
err = norm(u-exp(lambda*x),inf);    % measure the error
FS = 'fontsize';
xlabel('x',FS,12)
ylabel('exp(x)',FS,12)
title(sprintf('Solution of IVP for exp(x) -- error = %7.2e',err),FS,14) 

%%
% To solve IVPs, Chebfun uses a time-stepping method rather than a global 
% spectral method. If we add a `cheobppref` object and turn the display on
options = cheboppref();
options.display = 'iter';

%%
% we get the following message when solving the IVP
u = mldivide(L, 0, options);

%%
% Chebfun identifies that the boundary conditions correspont to an IVP and
% automatically uses a time-stepping method.
