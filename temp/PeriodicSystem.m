%% A periodic ODE system
% Nick Hale, October 2014

%%
% (Chebfun example ode-linear/PeriodicSystem.m)
% [Tags: #linearODE, #ODEsystem, #periodic, #piecewise]

%%
% Chebfun can solve systems of ODEs with periodic boundary conditions.
% For example, consider the equations
%
% $$ u - v' = 0, \qquad u'' + v = \cos(x) $$
%
% on the interval $[-\pi, \pi]$ with periodic boundary conditions on $u$ and
% $v$. A Chebfun solution could be put together like this:
d = [-pi,pi];
A = chebop(d);
A.op = @(x,u,v) [u-diff(v); diff(u,2)+v];
x = chebfun('x',d);
f = [0;cos(x)];
A.bc = 'periodic';
u = A\f;
u{1}, u{2}

%%
% Because the boundardy conditions are periodic, the system of ODEs is solved 
% with a Fourier collocation method, and the solution $u$ is represented by
% a Fourier series. (This is what `periodic` means in the display of $u$
% above.) We plot the result:
LW = 'linewidth'; lw = 2; FS = 'fontsize'; fs = 14;
plot(u,LW,lw), title('Solutions u and v',FS,fs), legend('u','v');

%%
% For this problem, the solution can actually be computed analytically.
% How close were we?
exact = [cos(x+3*pi/4) cos(x+pi/4)]/sqrt(2);
err = max([norm(u{1}-exact(:,1),inf) norm(u{2}-exact(:,2),inf)])

%%
% We show this also works for piecewise problems by artificially
% introducing a breakpoint at the origin.

A.domain = [-pi,0,pi];
u = A\f;
u{1}, u{2}

%%
% The solution is now represented by a Chebyshev series, and the equation
% has been solved with a Chebyshev collocation method, because Fourier
% collocation methods don't work with breakpoints.
plot(u,LW,lw), title('Solutions u and v',FS,fs), legend('u','v');
err = max([norm(u{1}-exact(:,1),inf) norm(u{2}-exact(:,2),inf)])
