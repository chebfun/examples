%% Allen-Cahn metastability
% Nick Trefethen, 25th November 2013

%% 
% (Chebfun Example pde/AllenCahn2.m)
% [Tags: #Allen-Cahn, #metastability, #EQUI, #PDE15S]

%% 1. The metastability phenomenon
% Oxford's Numerical Analysis Problem Solving Squad faced a problem this month involving
% the 1D Allen-Cahn equation, which is the following nonlinear PDE: $$ u_t =
% u_{xx} + u(1-u^2) . $$ This equation has three constant solutions: $u(x)=-1$,
% $u(x)=1$, and $u(x) = 0$, the first two stable and the third unstable.
% Depending on initial and boundary conditions, solutions tend to develop
% structures that locally approximate the value $1$ or $-1$, but as time passes,
% these structures eventually simplify. One of the built-in demos in CHEBGUI
% illustrates this behavior.

%%
% For example, if $x$ ranges over the whole real line with boundary conditions
% $u(x) \to -1$ as $|x|\to \infty$, here is the solution for the initial
% function $$ u(x) = 3 e^{-x^2/4} - 1, $$ computed with Chebfun's PDE15S
% command. Thanks to the exponential decay, we can approximate the real line by
% the finite interval $[-16, 16]$.
u0 = chebfun('3.*exp(-x.^2/4)-1',[-16 16]);
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
opts = pdeset('Eps',1e-6,'plot','off');
pdefun = @(u, t, x, diff) diff(u, 2)+u-u.^3;
t50 = 0:50; bc.left = -1; bc.right = -1;
[t, u] = pde15s(pdefun, t50, u0, bc, opts);
clf, waterfall(u, t, 'simple', LW, 1), grid on, 
view([115 20]), colormap([0 0 0]), xlabel('t', FS, 14)

%%
% Look what happens if we double the width of the initial pulse:
u0 = chebfun('3.*exp(-x.^2/16)-1', [-16 16]);
[t u] = pde15s(pdefun, t50, u0, bc, opts);
waterfall(u, t, 'simple', LW, 1), grid on, view([115 20]), colormap([0 0 0])
xlabel('t', FS, 14)

%%
% The time scale has become MUCH longer as we can see by increasing the time
% axis by a factor of 200:
u0 = chebfun('3.*exp(-x.^2/16)-1', [-16 16]);
t10000 = 0:200:10000;
[t, u] = pde15s(pdefun, t10000, u0, bc, opts);
waterfall(u, t, 'simple', LW, 1), grid on, view([115 20]), colormap([0 0 0])
xlabel('t', FS, 14)

%% 2. Computing a critical time
% The Problem Squad was asked, for the initial condition
u0 = chebfun('3.*exp(-x.^2/6)-1', [-16 16]);

%%
% what is the time at which the solution becomes negative?  A quick run suggests
% that the answer will be a little greater than 50:
t100 = 0:2:100;
[t, u] = pde15s(pdefun, t100, u0, bc, opts);
waterfall(u, t, 'simple', LW, 1), grid on, view([115 20]), colormap([0 0 0])
xlabel('t', FS, 14)


%%
% Let's see how the maximum decays with time:
umax = max(u(:, :));
plot(t100, umax, '.', MS, 14), grid on
xlabel('t', FS, 14), ylabel('max(u)', FS, 14) 

%%
% Where does the curve pass through zero?  All kinds of interpolation methods
% could be used to get an approximate answer to that question.  For fun, let us
% globally interpolate these equispaced data using the 'equi' flag in the
% Chebfun constructor:
p = chebfun(umax, [min(t) max(t)], 'equi');
plot(p, LW, 1.6), grid on
xlabel('t', FS, 14), ylabel('max(u)', FS, 14) 

%%
% The Gibbs wiggles show that for this coarsely sampled data set, this is not an
% excellent interpolant. Nevertheless, let's use it to estimate the critical
% time:
critical_time = roots(p)

%%
% This is actually pretty good, matching to several digits the estimate from a
% more careful computation, $55.92792$.
