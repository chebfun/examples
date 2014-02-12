%% Time-dependent integro-differential equation
% Nick Hale, October 2010

%%
% (Chebfun example pde/IntegroDiffT.m)
% [Tags: #integrodifferential]

%%
% Here we demonstrate how to solve the time-dependent integro-differential
% equation 
%
%  u_t = 0.02 u"(x) + \int_{-1}^{1}u(\xi)d\xi \int_{-1}^{x}u(\xi)d\xi,
%  u(-1) = u(1) = 0.
%
% using Chebfun's PDE15S command.

%%
% We work on the domain [-1,1] with a pulse initial condition:
d = [-1,1];
x = chebfun('x',d);
u0 = (1-x.^2).*exp(-30*(x+.5).^2);

%%
% Here is the anonymous function defining the problem for PDE15S.
f = @(u,t,x,diff,sum,cumsum) 0.02*diff(u,2) + cumsum(u)*sum(u);

%%
% The 4th, 5th, and 6th arguments define the differential, integral 
% (sum), and indefinite integral (cumsum) operators, respectively. 
% See 'help pde15s' for more details.

%%
% Now we solve the problem and plot the result.
[tt uu] = pde15s(f,0:.1:4,u0,'dirichlet');
FS = 'fontsize';
waterfall(uu,tt,'simple')
xlabel('x',FS,14)
ylabel('t',FS,14)
zlabel('z',FS,14)

%%
% This example can also be found as the "Integro-differential equation"
% demo among the PDE-Scalar demos of CHEBGUI.
