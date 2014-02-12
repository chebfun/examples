%% AN ODE ON AN UNBOUNDED INTERVAL
% Nick Hale, November 2010

%%
% (Chebfun example ode/UnboundedODE.m)

%%
% Chebfun has some support for ODEs on unbounded intervals (although this 
% feature is still in the experimental stages!).
%
% Here we solve equation 
%
%    0.1 u'' + u' + u = 0 ,    u(0) = 1,   u(inf) = 0;

[d x] = domain(0,inf);
A = -0.1*diff(d,2) + diff(d) + eye(d);
A.lbc = 1;
A.rbc = 0;
u = A\0;

LW = 'linewidth'; lw = 2; FS = 'fontsize'; fs = 18;
plot(u,LW,lw)
title('Solution',FS,fs)

%%
% We can check the solution by plotting the residual
% close, plot(A*u,LW,lw)
% title('Residual',FS,fs)

%% 
% We can also solve piecewise problems on unbounded domains. As an example,
% we'll take the same equation, and artificial insert a breakpoint at x = 1

[d x] = domain(0,inf);
B = -0.1*diff(d,2) + diff(d) + eye(d) - 1.5*diag(x<2);
B.lbc = 1;
B.rbc = 0;
v = B\0;

%%
% The solution should of course be the same:
figure
plot(v,'r',LW,lw)
title('Solution',FS,fs)
norm(u-v,inf)
