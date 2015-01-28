%% Nearest orthonormal functions
% Behnam Hashemi, December 2014

%%
% (Chebfun Example approx/NearestOrthFun.m)
% [Tags: #orthogonal polynomials, #quasimatrix, #SVD, #gallery]

function nearestOrthFun()
%%
% Suppose $A$ is a given $m \times n$ matrix. The problem of finding the
% orthonormal matrix $Q$ nearest to $A$ is well-known [4] (By an orthonormal
% matrix, we mean a matrix with orthonormal columns.). Specifically,
% it is the following problem:
%    $$ \mbox{argmin}\_U ||A-U||_\mbox{fro} \quad\mbox{subject to}\quad U^T U = I. $$
% This problem is a special case of the orthogonal Procrustes problem [3, p. 327].
% There are different ways to find the unique solution to this problem.
% One way that expresses $Q$ explicitly requires finding a matrix inverse
% square root [5] as follows:
%
% $$ Q = A(A^T A)^{-1/2}. $$
%
% However, the simplest (and possibly a numerically more stable) way is to
% compute the singular value decomposition of $A$ and to replace the singular
% values with ones. In other words, if $A = U S V^T$, then
% $Q = U V^T$ is the nearest orthonormal matrix to $A$. This is
% actually the unitary factor in the polar decomposition of $A$ [4].
%
% Our goal is to generalize this situation to the continuous case.
% Suppose that a set of functions defined on a domain are given.
% We know how to use Chebfun's `qr` command to compute orthogonal
% functions (stored as the columns of a qusimatrix `Q2`) corresponding
% to our given functions. See e.g. [2, Ch. 6] for more details. Here, we
% are looking for the nearest set of orthonormal functions. In Chebfun
% this is also straightforward. We just
% need to create a quasimatrix out of our given set of functions and use
% Chebfun's overloaded `svd` command [1]. The function
% `[Q,Q2] = ortho(A)` in the following computes both `Q` and `Q2`.

%% The first six monomials
A = cheb.gallery('vandermonde'); %A = [1, x, x.^2, x.^3, x.^4, x.^5];
[Q,Q2] = ortho(A);

%%
% Here, we see Legendre polynomials scaled to have norm 1 (left)
% and the nearest set of orthonormal functions (right). The first
% polynomial among these nearest orthonormal
% polynomials is approximately the following: $ p(x) = 0.9200 - 0.9035 x^2 +
% 0.3084 x^4$.

%% Chebyshev-Vandermonde quasimatrix on the interval [0, 1]
A = cheb.gallery('vandercheb');
A = restrict(A,[0,1]);
[Q,Q2] = ortho(A);
%% Another example
x = chebfun('x'); A = [1, cos(x), sin(x.^2), x.^3, x.^4, x.^5];
[Q,Q2] = ortho(A);
%% Non-smooth functions on the interval [0, 10]
A = [cheb.gallery('stegosaurus'), cheb.gallery('wiggly'), cheb.gallery('blasius')];
[Q,Q2] = ortho(A);
%%
function [Q,Q2] = ortho(A)
[U,S,V] = svd(A);
Q = U*V'; % nearest quasimatrix to A with orthonormal columns

[Q2,R] = qr(A);

m = size(A,2);
figure;
subplot(1,2,1);
plot(Q2); grid on; title('QR orthonormalization')
subplot(1,2,2);
plot(Q); grid on; title('Optimal orthonormalization')
fprintf(['Departure from orthogonality in the columns' ...
    ' of A = %3.2f\n'], norm(A'*A-eye(m)))
fprintf(['Departure from orthogonality in the columns' ...
    ' of Q = %3.2g\n'], norm(Q'*Q-eye(m)))
fprintf(['Departure from orthogonality in the columns' ...
    ' of Q2 = %3.2g\n'], norm(Q2'*Q2-eye(m)))
fprintf('The distance between A and Q2 = %3.2f\n', norm(A-Q2,'fro'))
fprintf(['The distance between A and its closest orthonormal' ...
    ' quasimatrix = %3.2f\n'], norm(A-Q,'fro'))
end

%% References
%
% 1. Z. Battles, L.N. Trefethen, An extension of MATLAB to continuous
%    functions and operators, SIAM J. Scientific Computing 25 (2004)
%    1743-1770.
%
% 2. T.A. Driscoll, N. Hale, and L.N. Trefethen, editors, _Chebfun Guide_,
%    Pafnuty Publications, Oxford, 2014.
%
% 3. G.H. Golub, C.F. Van Loan, Matrix Computations, 4th edition, The Johns
%    Hopkins University Press, 2013.
%
% 4. N.J. Higham, Matrix nearness problems and applications, pp. 1-27, in
%    Applications of Matrix Theory, Inst. Math. Appl. Conf. Ser. New Ser.,
%    Vol. 22, Oxford Univ. Press, New York, 1989.
%
% 5. B.K.P. Horn, H.M. Hilden, S. Negahdaripour, Closed-form solution of
%    absolute orientation using orthonormal matrices, Journal of the Optical
%    Society A, 5 (1988) 1127-1135.
end