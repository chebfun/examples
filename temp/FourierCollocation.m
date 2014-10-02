%% Fourier Spectral Collocation 
% Hadrien Montanelli, October 2014

%%
% (Chebfun example temp/FourierCollocation)
% [Tags: #linearODE, #periodic]
LW = 'linewidth'; dom = [0 2*pi];

%%
% A Fourier spectral collocation method is now available in Chebfun to
% solve ODEs with periodic boundary conditions. The solution is a chebfun
% using a `fourtech` representation, that is, a trigonometric interpolant on 
% equispaced points. 
% This is the default method for periodic boundary conditions.
%
% Consider the following first-order ODE
%
% $$ u'(x) + a(x)u(x) = f(x) $$
%
% on $[0,2\pi]$, with periodic boundary conditions, and where $a(x)$ and $f(x)$  
% are continuous and periodic complex-valued functions. 
% This equation has a unique periodic solution if 
% $\overline{a}=\frac{1}{2\pi}\int_0^{2\pi}a(x)dx\neq ik$ for all integers k. 
% In particular, if $a(x)=a$ is a constant coefficient, this means 
% $a\neq ik$ for all $k$. 
%
% Take for example $a(x)=1+\sin(\cos(10x))$ and $f(x)=\exp(\sin(x))$, and solve 
% it with Fourier collocation. Since $\overline{a}=1$, this a well-posed
% problem.
L = chebop(@(x,u) diff(u) + (1+sin(cos(10*x))).*u, dom); 
L.bc = 'periodic';
f = chebfun(@(x) exp(sin(x)), dom);
u = L \ f
figure, plot(u, LW, 2)

%%
% The periodic solution $u$ satisfies the differential equation to high
% accuracy:
norm(L*u - f, inf)

%% 
% We can solve the same ODE with Chebyshev collocation on 2nd-kind points using 
% a `cheboppref` object with `chebcolloc2` discretization. (Chebyshev
% collocation on 1st-kind points is also possible, use `chebcolloc1`.)
pref = cheboppref();
pref.discretization = 'chebcolloc2';
v = solvebvp(L, f, pref)
hold on, plot(v, 'r', LW, 2)

%%
% The solution $v$ is now a chebfun with a `chebtech2` representation, that is, 
% a polynomial interpolant on on 2nd-kind Chebyshev points.
% It satisfies the differential to high accuracy too
norm(L*v - f, inf)

%%
% but is about $\pi/2$ times longer.
length(v)/length(u)

%%
% Consider now the second-order ODE 
%
% $$ u''(x) +  a_1(x)u'(x) + a_0(x)u(x) = f(x) $$
%
% on $[0,2\pi]$, with periodic boundary conditions, and where $a_0(x)$,
% $a_1(x)$, and $f(x)$ are continuous and periodic complex-valued functions.
% Let $\Delta$ be the Hill discriminant of this equation
%
% $$ \Delta = \frac{c(2\pi) + s'(2\pi)}{2}, $$
%
% where $c(x)$ and $s(x)$ are the solutions of the homogeneous equation, 
% corresponding to the initial conditions $c(0)=1$, $c'(0)=0$ and $s(0)=0$, 
% $s'(0)=1$. The nonhomogeneous equation has a unique periodic solution if 
% $\Delta \neq 1$ [1].
%
% Take $a_1(x)=\sin(\cos(x/2)^2)$, $a_0(x)=\cos(12\sin(x))$, and 
% $f(x)=\exp(\cos(2x))$, and solve it with Fourier collocation.
a1 = chebfun(@(x) sin(cos(x/2).^2), dom);
a0 = chebfun(@(x) cos(12*sin(x)), dom);
L = chebop(@(u) diff(u, 2) + a1.*diff(u) + a0.*u, dom); 
L.bc = 'periodic';
f = chebfun(@(x) exp(cos(2*x)), dom);
u = L \ f
figure, plot(u, LW, 2)

%%
% Again, the periodic solution $u$ satisfies the differential equation to high 
% accuracy
norm(L*u - f, inf)

%%
% The solution with Chebyshev collocation on 2nd-kind points 
pref = cheboppref();
pref.discretization = 'chebcolloc2';
v = solvebvp(L, f, pref)
hold on, plot(v, 'r', LW, 2)

%%
% is about 2.35 times longer:
length(v)/length(u)

%% 
% The second-order ODE we have solved is well-posed, and we can check that 
% computing the Hill discriminant, and verifying that it is not 1:
L.bc = [];
L.lbc = @(c) [ c - 1 ; diff(c) ];
c = L \ 0;
L.lbc = @(s) [ s ; diff(s) - 1 ];
s = L \ 0;
HillDiscr = 1/2*(c(2*pi) + feval(diff(s), 2*pi))

%%
% If we take $a_1(x) = \sin(x)$ and $a_0(x)=1+ia_1(x)$, the solution we
% obtain does not satisfy the differential equation:
L = chebop(@(x,u) diff(u, 2) + sin(x).*diff(u) + (1+1i*sin(x)).*u, dom);
L.bc = 'periodic';
u = L \ f;
norm(L*u - f, inf)

%%
% This correponds to a Hill discriminant equal to 1 (up to 14 digits):
L.bc = [];
L.lbc = @(c) [ c - 1 ; diff(c) ];
c = L \ 0;
L.lbc = @(s) [ s ; diff(s) - 1 ];
s = L \ 0;
HillDiscr = 1/2*(c(2*pi) + feval(diff(s), 2*pi))

%% References
%
% 1. M. S. P. Eastham, _The spectral theory of periodic differential
%    equations_, Scottish Academic Press, 1973.