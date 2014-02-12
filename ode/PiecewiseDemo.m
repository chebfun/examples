%% Piecewise operators demo
% Nick Hale, 25th November 2011

%%
% (Chebfun Example ode/PiecewiseDemo.m)
% [Tags: #linearODE, #piecewise]

%%
% Here we demonstrate piecewise differential operators (incl. boundary
% conditions), and how the systems constructor goes about solving them in
% Chebfun. This Example is intended for those wanting to learn more about
% what's going on under the cheb-hood, rather than those simply wanting to
% solve these kinds of problems. If you fall into this later category,
% there are more Examples of that kind in the ODE examples directory [1].

format short

%% Piecewise operator

%%
% Define a chebop with coefficient which jumps at x = 0.

A = chebop(@(x,u) -diff(u,2) + sign(x).*exp(x).*u)

%%
% Although A doesn't know it has a discontinuous coefficient at zero, it
% will once we perform a linearity check (which is done internally when we
% solve a BVP).

A = linop(A)

%%
% This would still be the case if A were nonlinear and we were linearising
% around a current iteration.

%%
% What does it look like when we evaluate these piecewise operators? Well
% without boundary conditions, with just have two independant blocks, which
% can be evaluated at different sizes (here only  4x4 for x in [-1 0] and
% 3x3 in [0 1], so that it fits on the screen!).

A([4 3],'nobc')

%%
% By default, Chebfun will apply boundary conditions to enforce continuity
% of derivatives up to the differential order of the operator. This can be
% seen in the bottom two rows below.

A([4 3],'bc')

%%
% However, we still need to apply some boundary conditions of our own to 
% the operator. Let's choose dirichlet for simplicity.

B = A & 'dirichlet';
B([4 3],'bc')

%%
% We'll need a RHS to solve for. Again, for simplicity let's choose the
% constant function 1. With the boundary conditions tagged on, for a given
% N this is then be given by

rhs = @(N) [ones(N-4,1) ; zeros(4,1)];

%%
% We're now almost in a position to start solving piecewise ODEs. However,
% the standard constructor doesn't do quite enough here, as when it
% constructs on a domain such as [-1 0 1], the 2 subdomains are treated
% independently. By wrapping the domain as a cell, we force the use of the
% systems constructor which doesn't suffer from this.

myfun = @(x,N,bks) B(N{:},'bc')\rhs(sum(N{:}));
u = chebfun(myfun,{[-1 0 1]},'eps',1e-10)

%%
% An alternative notation is
u = chebfun(myfun,[-1 0 1],'eps',1e-10,'sys',1);

%%
% Of course most users won't even see things at this level - they'll just 
% be calling backslash!
v = B\1;

%%
% which we see does much the same as above.
plot(u,'b',v,'--r','LineWidth',1.6)

%% Piecewise RHS
% Jumps and discontinuities can also be introduced by the RHS. For example,
% here's a smooth operator

A = chebop(@(x,u) -diff(u,2) + exp(x).*u);
x = chebfun('x');
A.bc = 0

%%
% and a RHS with a jump at sqrt(2)/2

s = sqrt(2)/2;
rhs = heaviside(x-s);

%%
% Again, when we solve the system, continuiuty of the solution and its
% first derivative are enforced across the discontinuity of the RHS at
% sqrt(2)/2.

v2 = A\rhs
plot(v2)

%% References
% [1] http://wwww.chebfun.org/examples/ode/
%
% [2] http://wwww.chebfun.org/examples/ode/html/JumpConditions.shtml



