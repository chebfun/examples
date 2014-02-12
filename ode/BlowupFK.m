%% Blowup equation (Frank-Kamenetskii)
% Nick Trefethen, 25 September 2010

%%
% (Chebfun example ode/BlowupFK.m)
% [Tags: #nonlinearODE, #blowup]

%%
% The Frank-Kamenetskii or "spontaneous combustion" equation is
% the PDE du/dt = d2u/dx2 + Aexp(u).  On the interval [-1,1] with
% zero initial and boundary conditions, solutions to this
% equation blow up to infinity in finite time if 
% A is bigger than about 0.878.  For smaller A, solutions
% converge to a steady state.

%%
% Here we compute some of these steady-state solutions, which
% are solutions of the ODE boundary value problem
% u"+A*exp(u)=0, u(-1)=u(1)=0.

N = chebop([-1 1]);
N.bc = 'dirichlet';
FS = 'fontsize'; figure
for A = [.2 .4 .6 .8 .87]
  N.op = @(u) diff(u,2) + A*exp(u);
  u = N\0;
  plot(u,'linewidth',2), grid on, hold on
  text(-.1,max(u)+.04,['A = ' num2str(A)],FS,14)
end
axis([-1 1 0 1.2])
title('Frank-Kamenetskii blowup equation',FS,16)

%%
% Reference:
%
% H. Fujita, On the nonlinear equations Del u + exp(u) = 0
% and dv/dt = Del v + exp(v), Bulletin of the American Mathematical
% Society 75 (1969), 132-135.
