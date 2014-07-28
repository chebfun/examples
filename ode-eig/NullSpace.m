%% The nullspace of a linear operator
% Nick Hale and Stefan Guettel, December 2011

%%
% (Chebfun example ode-eig/NullSpace.m)
% [Tags: #nullspace, #NULL, #SUBSPACE]


%%
% We've recently introduced some new functionality in Chebfun for computing the
% nullspace of differential operators. Let's explore this with a couple of
% simple examples.

%% 1. Simple example #1
% Let's start as simply as we can, and take
%
% $$ (Lu)(x) := u''(x), \quad x\in [-1, 1]. $$
L = chebop(@(u) diff(u, 2));

%%
% Clearly the nullspace of this operator -- that is, the space
% of functions $v$ for which
% $L(v)=0$ -- is spanned the two functions
v = [1, chebfun('x')];
norm(L(v))

%%
% Supposing we didn't know this, we could compute a basis
% for the nullspace with the `null` method:
LW = 'LineWidth'; lw = 1.6;
V = null(L)
plot(V, LW, lw)
V'*V
norm(L(V))

%%
% where we find that $V^T V = I$ and $LV \approx 0$ as required.

%%
% Clearly `V` doesn't correspond directly to $1$ and $x$, since there is some
% freedom in how we orthogonalise the basis. However, we can check that `V` and
% $\{1, x\}$ correspond to the same spaces by computing the angle between the
% spaces with the `subspace` command.

subspace(v, V)

%% 2. Incomplete boundary conditions
% Now let's consider the more complicated 2nd-order operator
%
% $$ Lu = u'' + 0.1x(1-x^2)u' + \sin(x)u, \quad
% x\in [-\pi, \pi].  \qquad (*) $$
dom = [-pi, pi];
L = chebop(@(x, u) diff(u, 2) + .1*x.*(1-x.^2).*diff(u) + sin(x).*u, dom);

%%
% As before, it has a nullspace of rank 2.
V = null(L)
plot(V, LW, lw)
V'*V
norm(L(V))

%%
% However, now suppose we impose one boundary condition, say,
% a Dirichlet condition at the
% left. This removes one degree of freedom, and we are left with a rank 1
% nullspace.
L.lbc = 0;
L.rbc = [];
v = null(L)
plot(v, LW, lw), shg
v'*v
norm(L(v))

%%
% Clearly this null vector must satisfy the given condition $v(-\pi) = 0.$
v(-pi)

%% 3. An application
% Where might these ideas be useful?
% Well, suppose we were interested in equation $(*)$
% with a homogeneous Dirichlet condition at the left, and wanted to know what
% inhomogeneous Dirichlet condition gave the minimal 2-norm of the solution to
% $Lu = 1$.
% Rather than solving the linear system for a number of different boundary
% conditions (which would be computationally expensive) we could simply solve
% for one, say again a homogeneous Dirichlet condition,
L.rbc = 0;
u = L\1;
hold on, plot(u, '--r', LW, lw), hold off

%%
% and compute the rest by adding a scalar multiple of the null-function $v$.
E = chebfun(@(c) norm(u + c*v, 2), [-10, 10], 'vectorize', 'splitting', 'on');
plot(E,LW,lw)

%%
% We compute the 2-norm as a chebfun in the unknown variable $c$, which we can
% then minimise to obtain the minimal energy solution
[minE, c_star] = min(E)
u_star = u + c_star*v
plot(u_star,LW,lw)

%%
% So the condition we require is that $u(\pi)$ = `bc_star`, where

bc_star = u_star(pi)
%% 4. Exotic constraints
% The Chebfun `null` function can also
% handle the more exotic types of boundary conditions
% that can be imposed in Chebfun (see [1]). For example, suppose we wish to
% compute the nullspace of the 3rd-order piecewise-smoooth ODE
%
% $$ Lu := 0.1u''' + \sin(x)u'' + u, \quad x\in[-1,1] $$
%
% with the 'boundary' condition
%
% $$ \int(u) = u(0). $$
dom = [-1, 1];
L = chebop(@(x, u) .1*diff(u, 3) + sin(x).*diff(u, 2) + u, dom);
L.bc = @(x, u) sum(u) - u(0);

%%
% Here `null` has no problems!
V = null(L)
plot(V, LW, lw), shg
V'*V

%%
sum(V) - V(0,:)
norm(L(V), 1)

%% 5. References
%
% 1. Chebfun Example [ode-linear/NonstandardBCs](../ode-linear/NonstandardBCs.html)
