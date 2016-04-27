%% Gray-Scott equation in 2D
% Nick Trefethen, April 2016

%%
% (Chebfun Example pde/GrayScott.m)
% [Tags: #Gray-Scott, #spin2]

%% 1. Spots
% The Gray-Scott equations are a pair of coupled reaction-diffusion
% equations that lead to interesting patterns.
% Let us look at some examples in 2D.

%%
% The equations are
% $$ u_t = \varepsilon_1\Delta u + b(1-u) - uv^2, \quad
% v_t = \varepsilon_2\Delta v - dv + uv^2, $$
% where $\Delta$ is the Laplacian and $\varepsilon_j,b,d$ are parameters.
% To begin with we choose these parameters:
ep1 = 0.00002; ep2 = 0.00001;
b = 0.035; d = 0.095;
%%
% We now solve up to $t=1000$.
tic, dom = [-1 1 -1 1]; x = chebfun('x'); tspan = [0 1000];
S = spinop2(dom,tspan);
S.linearPart = @(u,v) [ep1*lap(u); ep2*lap(v)];
S.nonlinearPart = @(u,v) [b*(1-u)-u.*v.^2;-d*v+u.*v.^2];
S.init = chebfun2v(@(x,y) 1-exp(-60*(x.^2+y.^2)), ...
                   @(x,y) exp(-60*(x.^2+(y+.01).^2)));
pause off
tic, u = spin2(S,spinpref2('plot','off')); t = toc;
plot(u{1}), view(0,90), axis equal
%FS = 'fontsize'; text(42,3.4,'t=0 and t=100',FS,12)

%% 3. References
%
% [1] _PDE Coffee Table Book_.
