%% Piecewise operators demo
% Nick Hale, November 2011
% Toby Driscoll, 9 June 2014

%%
% (Chebfun Example ode/PiecewiseDemo.m)
% [Tags: #linearODE, #piecewise]

%%
% Here we demonstrate how piecewise differential operators are discretized
% in Chebfun. These details are not necessary to know in order to use
% Chebfun to solve such problems. Rather, they are to satisfy the curious
% about what is happening under the hood. 

%%
format short
plotopt = {'linewidth',2,'markersize',10};
cheboppref.setDefaults('discretization',@colloc2);

%% Piecewise operator

%%
% Define a chebop with coefficient which jumps at x = 0.

L = chebop(@(x,u) -diff(u,2) + sign(x).*exp(x).*u);

%%
% We'll add homogenous Dirichlet conditions to the operator.
L = L & 'dirichlet';

%%
% Although A doesn't know it has a discontinuous coefficient at zero, it
% will once we perform a linearity check (which is done internally when we
% solve a BVP).

A = linop(L);   A.domain

%%
% This would still be the case if A were nonlinear and we were linearising
% around a current iteration.

%%
% What are appropriate discretizations of such a piecewise operator? 
% To begin with, let's see what happens when we have just 2 points in [-1,0]
% and 3 in [0,1].

matrix(A,[2 3])

%%
% The first two rows show the boundary conditions being applied at each
% endpoint. The remainder of the matrix consists of a 2x4 block and a 3x5
% block. These are rectangular because the operator is of second order, and
% a second derivative evaluated at 3 points (polynomial of degree 2) should
% be applies to a function evaluated at 5 points (polynomial of degree 4). 

%%
% The key observation, though, is that the first 4 columns are completely
% decoupled from the last 5. This can't give the right answer to the
% original, global problem! The appropriate coupling can be derived
% automatically by Chebfun.
B = deriveContinuity(A);
matrix(B,[2 3])

%%
% Now the first two rows contain the discretizations of the conditions 
% $u(0^+) - u(0^-) = 0$, $u'(0^+) - u'(0^-) = 0$. They couple together the
% otherwise separate blocks:
spy(B)

%%
% The two continuity rows, together with the two boundary condition rows,
% make up for the rectangularity of the block operators, so that the entire
% matrix is square. 
size( matrix(B,[2 3]) )

%% Discrete problem

%%
% We can now see how to solve a discrete version of the problem L(u)=1. We
% have to discretize the function 1 on n points in each subinterval, and
% prepend four zeros for the boundary and continuity conditions.
n = 12;
b = [ zeros(4,1); ones(2*n,1) ];

%%
% Then we solve the discrete system. 
u = matrix(B,[n n]) \ b;

%%
% The solution is defined on two 'copies' of the Chebyshev points.
% Remember, the solution is defined on two extra points in each subinterval
% due to the rectangularity of differentiation.
x0 = chebtech2.chebpts(n+2);
x = [ (x0-1)/2; (x0+1)/2 ];
plot(x, u, 'o-', plotopt{:} )
hold on
plot(x, 0*u-0.05, '.k' , plotopt{:})
grid on

%% Automatic solution

%%
% All of the continuity manipulations happen behind the scenes when you
% solve the BVP using backslash (or an eigenvalue problem with EIGS).
u = L\1
plot(u,'r',plotopt{:})

%%
% We can verify the continuity conditions:
j0 = jump(u,0)
j1 = jump(diff(u),0)

% Jumps can also be introduced by the right-hand side function of the BVP.
% For example, here's a smooth operator.

A = chebop(@(x,u) -diff(u,2) + exp(x).*u) & 'dirichlet';

%%
% and a RHS with a jump at sqrt(2)/2

x = chebfun('x');
f = heaviside( x - 1/sqrt(2) );
clf, plot(f*.021,'k',plotopt{:})

%%
% Again, when we solve the system, continuiuty of the solution and its
% first derivative are enforced across the discontinuity of the RHS at
% 1/sqrt(2).

v = A\f

%%
hold on
plot(v,plotopt{:})

%%
% Because the function on the right side has a jump, and the solution is
% continuous, the second derivative of the solution must also have a jump.
clf, plot((f-0.5)*.41,'k',plotopt{:})
hold on, plot(diff(v,2),'r',plotopt{:})
