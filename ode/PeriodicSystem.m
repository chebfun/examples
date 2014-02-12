%% A periodic ODE system
% Nick Hale, November 2010

%%
% (Chebfun example ode/PeriodicSystem.m)
% [Tags: #linearODE, #ODEsystem, #periodic, #piecewise]

%%
% Chebfun can solve systems of ODEs with periodic boundary conditions.
% For example, consider the equation 
%
%   u  -  v' = 0,
%   u" +  v  = cos(x),
%
% on the interval [-pi, pi] with periodic boundary conditions on u and v.
% A Chebfun solution could be put together like this:

d = [-pi,pi];
A = chebop(d);
A.op = @(x,u,v) [u-diff(v), diff(u,2)+v];
x = chebfun('x',d);
f = [0, cos(x)];
A.bc = 'periodic';
u = A\f;

%%
% We plot the result:
LW = 'linewidth'; lw = 2; FS = 'fontsize'; fs = 14;
plot(u,LW,lw), title('Solutions u and v',FS,fs), legend('u','v');

%%
% For this problem, the solution can actually be computed
% analytically.  How close were we?
true = [cos(x+3*pi/4) cos(x+pi/4)]/sqrt(2);
err = norm(u-true,inf)

%%
% We show this also works for piecewise problems by artificially
% introducing a breakpoint at the origin.

A.domain = [-pi,0,pi];
u = A\f;
plot(u,LW,lw), title('Solutions u and v',FS,fs), legend('u','v');
err = norm(u-true,inf)
