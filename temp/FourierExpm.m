%% Time-dependent PDE on a periodic interval
% Hadrien Montanelli, October 2014

%%
% (Chebfun example temp/FourierEigs)
% [Tags: #heatequation, #convection, #periodic, #linearPDE, #EXPM]
LW = 'linewidth'; dom = [0 2*pi];

%%
% Consider the time-dependent PDE on $[0,2\pi]\times [0,T]$
%
% $$ u_t = \mathcal{L}u, $$
%
% on $[0,2\pi]\times[0,\infty)$, with given periodic initial condition 
% $u(x,0)$. We seek periodic solutions $u(x,t)$. If the operator 
% $\mathcal{L}$ is semi-bounded, this problem is well-posed. This is the 
% case for the heat equation, $\mathcal{L}u=u_{xx}$, and for the convection
% equation, $\mathcal{L}u=c(x)u_x$.

%%
% Consider first the heat equation
%
% $$ u_t = u_{xx}, $$
%
% on $[0,2\pi]\times[0, T]$, with periodic boundary conditions, and with
% initial condition $u(x,0)=\sin(3x)$.
% We can solve it in Chebfun as follows.
T = 1; dt = 0.05;
t = [0:dt:T];
L = chebop(@(u) diff(u, 2), dom);
L.bc = 'periodic';
u0 = chebfun(@(x) sin(3*x), dom); 
u = expm(L, t, u0);
figure, waterfall(u, t, LW, 2)
view(10, 70), axis([0 2*pi 0 T -1 1])

%%
% The diffusion has done its job: the solution at $T=1$ has very small
% amplitude.
norm(u{end}, 'inf')

%% 
% Let us solve now the convection equation
%
% $$ u_t = c(x)u_x, $$
%
% on $[0,2\pi]\times[0, T]$, with $c(x)= -\frac{1}{5}-\sin^2(x-1))$,
% periodic boundary conditions, and initial condition 
% $u(x,0)=\exp(-100(x-1)^2)$.
T = 8; dt = 0.5;
t = [0:dt:T];
c = chebfun(@(x) -(1/5 + sin(x-1).^2), dom);
L = chebop(@(x,u) c.*diff(u, 1), dom);
L.bc = 'periodic';
u0 = chebfun(@(x) exp(-100*(x-1).^2), dom); 
u = expm(L, t, u0);
figure, waterfall(u, t, LW, 2)
view(10, 70), axis([0 2*pi 0 T 0 1])
