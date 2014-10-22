%% Nonlinear Periodic ODE
% Hadrien Montanelli, October 2014

%%
% (Chebfun example temp/FourierCollocationNonLin)
% [Tags: #nonlinearODE, #periodic]
LW = 'linewidth'; dom = [0 2*pi];

%%
% Chebfun uses Fourier-collocation to solve linear, non linear and system 
% of ODEs with periodic boundary conditions. Consider the following 
% nonlinear ODE       
% 
% $$ u' - u\cos(u) = \cos(4x), $$
%
% on [0, 2*pi], with periodic boundary conditions. For nonlinear ODEs, we 
% need to specify an intial guess; let us try $\cos(x)$. 
% We can solve the ODE on Chebfun as follows. 
f = chebfun(@(x) cos(4*x), dom);
N = chebop(@(u) diff(u) - u.*cos(u), dom);
N.bc = 'periodic';
N.init = chebfun(@(x) cos(x), dom);
u = N \ f

%%
% Let us plot the initial guess in green, the right-hand side in red, and 
% the solution in blue.
figure, plot(N.init, 'g', LW, 2)
hold on, plot(f, 'r', LW, 2)
hold on, plot(u, 'b', LW, 2)

%%
% The solution $u(x)$ satisfies the ODE to high accuracy:
norm(N*u - f, inf)

%%
% If we start with an other initial guess, we obtain an other solution.
% Let us try $\sin(x)^2$, plot it in dashed green, and plot the solution
% in dashed blue.
N.init = chebfun(@(x) sin(x).^2, dom);
v = N \ f
hold on, plot(N.init, '--g', LW, 2)
hold on, plot(v, '--b', LW, 2)

%%
% The solution $v(x)$ satisfies the ODE to high accuracy too:
norm(N*v - f, inf)
