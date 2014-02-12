%% Spectral radius of the SOR iteration matrix
% Nick Trefethen, 16th October 2012

%%
% (Chebfun example linalg/SOR.m)

%%
% [Tags: #linearalgebra, #spectralradius, #SOR]

%%
% The classic finite-difference 1D Laplacian discretization looks like
% this:
N = 11;
A = toeplitz([2 -1 zeros(1,N-3)])

%%
% We may split $A$ into its lower-triangular, diagonal, and
% upper-triangular parts:
L = tril(A,-1);
D = diag(diag(A));
U = triu(A,1);

%%
% From the beginning of the computer era, people studied solution of matrix
% problems with this kind of matrix by the method of _successive
% overrelaxation_ or _SOR_.  Here $\omega\in [1,2]$ is the overrelaxation
% parameter, and we iterate with the matrix defined like this: $$ G =
% M^{-1} N, \qquad M = D + \omega L, \quad N = (1-\omega)D- \omega U. $$ In
% Matlab, that's
G = @(om) (D+om*L)\((1-om)*D-om*U);
rho = @(om) max(abs(eig(G(om))));

%%
% Analysis of the SOR iteration was carried out by Frankel [1] and
% generalized by Young [4]; see also [3].  Details are given in innumerable
% books, such as Golub and Van Loan [2]. Supposing we didn't know the
% theory, Chebfun would give us an elegant way to draw the famous
% optimal-omega curve:

%%
f = chebfun(rho,[1 2],'splitting','on','vectorize');
plot(f,'linewidth',2), grid on
xlabel('\omega','fontsize',16)
ylabel('convergence factor')

%%
% Chebfun gives us the following optimal omega:
[rho_opt,omega_opt] = min(f)

%%
% Here are the exact optimal values:
omega_exact = 2/(1+sin(pi/N))
rho_exact = omega_exact - 1

%%
% References:
%
% [1] S. Frankel, Convergence rates of iterative treatments of partial
% differential equations, Mathematics of Computation 4 (1950), 56-75.
% 
% [2] G. H. Golub and C. F. Van Loan, Matrix Computations, 3rd ed.,
% Johns Hopkins, 1996.
%
% [3] R. J. LeVeque and L. N. Trefethen, Fourier analysis of the SOR
% iteration, IMA Journal of Numerical Analyaisis 8 (1988), 273-279.
%
% [4] D. M. Young, Iterative Methodfs for Solving Partial Difference
% Equations of Elliptic Type, PhD thesis, Harvard U., 1950.