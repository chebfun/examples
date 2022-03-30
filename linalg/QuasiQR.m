%% Backward stability of quasimatrix QR factorization
% Nick Hale and Nick Trefethen, March 2022

%%
% (Chebfun example linalg/QuasiQR.m)
% [Tags: #QR]

%%
% Wilkinson invented backward error analysis to explain
% the behavior of matrix algorithms.  Here is a nice example of
% how his ideas apply also to the continuous analogue, namely
% quasimatrix algorithms.

%%
% We start from a quasimatrix whose columns are Runge functions with varying
% scale parameters.  A call to |plotcoeffs| shows that that each column
% is resolved by an exponentially decreasing Chebyshev series.
x = chebfun('x');
A = 1/(1 + (1:10).*(x-.1)^2);
plotcoeffs(A), title A

%%
% This quasimatrix is highly ill-conditioned, as we can see by computing
% its condition number:
cond(A)

%%
% What's particularly interesting is to compute the quasimatrix QR factorization
% of $A$ and to look at the Chebyshev coefficients of the columns of $Q$:
[Q,R] = qr(A);
plotcoeffs(Q), title Q

%%
% Rounding errors are affecting these results fundamentally.  The first
% column of $Q$ has a Chebyshev series that converges down to machine precision, but with each 
% successive column, about one digits is lost.  Mathematically, all these columns
% should be smooth, so clearly rounding errors have messed up the mathematical picture.
% The computed $Q$ probably matches the mathematical ideal to about 6 digits, not 16.

%%
% But now here comes Wilkinson.  Wilkinson's startling discovery was that even though $Q$ is
% full of errors, and so is $R$, the errors in these two objects are ``diabolically correlated'' so that
% when you take the product, it matches $A$ to full precision:
norm(A-Q*R)

%%
% As he liked to put it, though $Q$ and $R$ are far from the correct factors of $A$, they are
% the exactly (or almost exactly) correct factors of a slightly perturbed matrix.
