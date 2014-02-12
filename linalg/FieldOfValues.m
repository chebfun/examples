%% Field of values and numerical abscissa
% Nick Trefethen, November 2010

%%
% (Chebfun Example linalg/FieldOfValues.m)
% [Tags: #linearalgebra, #fieldofvalues]

%%
% If A is a matrix, the field of values F(A) is
% the nonempty bounded convex set in the complex plane
% consisting of all the Rayleigh quotients of A, that is,
% all the numbers q'Aq, where q is a unit vector and
% q' is its conjugate transpose.  

%%
% The standard method for computing the field of
% values numerically is an algorithm due to C. R.
% Johnson in 1978 based on finding the
% maximum and minimum eigenvalues of (A+A')/2, then
% "rotating" this computation around in the complex plane [1].
% This algorithm is implemented in the Chebfun command FOV,
% which is listed at the end of this Example.

%%
% Generically the boundary of the field of values is
% smooth, but it is not always smooth.  Chebfun's
% 'splitting' feature enables FOV to compute this boundary
% in either situation, smooth or not.

%%
% For example, here are the eigenvalues and
% field of values of a random matrix of dimension 20.
% This is a case where the boundary is smooth.
randn('seed',1), A = randn(20);
LW = 'linewidth'; lw = 1.6; MS = 'markersize'; ms = 18;
FA = fov(A);
figure, plot(FA,LW,lw), axis equal, grid on
ax = axis; axis(1.1*ax)
hold on, plot(eig(A),'.k',MS,ms)

%%
% The numerical abscissa of A is the maximum real part of
% its field of values:
[alpha,maxtheta] = max(real(FA))

%%
% Here we add it to the plot as a red dot:
plot(real(FA(maxtheta)),imag(FA(maxtheta)),'.r',MS,24)

%%
% You can also find the numerical abscissa without Chebfun:
alpha = max(eig((A+A')/2))

%%
% Now let's consider a matrix B defined as a
% diagonal matrix with the same eigenvalues as A.
% In this case the boundary of the field of values is a polygon:
B = diag(eig(A));
FB = fov(B);
hold off, plot(real(FB),imag(FB),'b',LW,lw,'jumpline',{'b',LW,lw})
hold on, plot(eig(B),'.k',MS,ms), axis(1.1*ax), axis equal, grid on
[alpha,maxtheta] = max(real(FB));
plot(real(FB(maxtheta)),imag(FB(maxtheta)),'.r',MS,24)

%%
% Since the field of values is not smooth, its boundary is
% a Chebfun with several pieces:
FB

%%
% Finally, here's an example where the boundary of the
% field of values mixes smooth curves with straight segments:
C = [0 3 0 0; -3 0 0 0; 0 0 0 3; 0 0 1 1]
FC = fov(C);
hold off, plot(real(FC),imag(FC),'b',LW,lw,'jumpline',{'b',LW,lw})
axis(4*[-1 1 -1 1]), axis square, grid on
hold on, plot(eig(C),'.k',MS,ms)

%%
% Here is a listing of FOV.
% Note that the numerical computations are carried out
% in just about 10 lines of code.
type fov

%%
% References:
%
% [1] C. R. Johnson, Numerical determination of the
% field of values of a general complex matrix,
% SIAM J. Numer. Anal. 15 (1978), 595-602. 
%
% [2] L. N. Trefethen and M. Embree, Spectra and Pseudospectra:
% The Behavior of Nonnormal Matrices and Operators, Princeton U.
% Press, 2005, chapter 17 on Numerical range, abscissa, and radius.

