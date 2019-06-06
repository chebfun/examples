%% Brownian bridge and random polynomials
% Nick Trefethen, June 2019

%%
% Chebfun example stats/BrownianPolynomials.m)
% [Tags: #randnfun, Brownian bridge]

%%
% The Chebfun |randnfun(d)| command produces a smooth random function
% with maximal wavelength $d$, as described in [1], 
% and its indefinite integral is
% a smooth random walk, which converges to Brownian motion
% as $d\to 0$.  If we subtract off the mean of the integrand, 
% we get a _smooth Brownian bridge_ starting and ending at 0, like this:
rng(1)
r = randnfun(0.02,[0,1]); r = r - mean(r);
b = cumsum(r);
plot(b), grid on

%%
% All this is done via trigonometric polynomials with random
% coefficients.  (The idea of Fourier series with random
% coefficients goes back to Norbert Wiener and is the subject
% of a marvelous book by Kahane [4].)  
% The resulting smooth random function distribution
% is effectively translation-invariant over its interval of definition.

%%
% A different approach to Brownian bridge has been introduced by James
% Foster and his coauthors
% recently in [2], based on ordinary algebraic polynomials with random
% coefficients.  The polynomials
% needed for the basis are like Jacobi polynomials,
% but with the nonstandard exponent $-1$ at both endpoints, and 
% constrained to take the value 0 there.  They
% can be computed for each $n\ge 2$ as a difference of
% two Legendre polynomials:
foster = @(n) ( legpoly(n,[0 1]) - legpoly(n-2,[0 1]) )/sqrt(8*n-4);

%%
% Here for example is the degree $50$ Foster polynomial:
plot(foster(50)), grid on

%%
% The same anonymous function works for multiple columns.  Here,
% for example, are the first 5 polynomials:
plot(foster(2:6)), grid on

%%
% The random polynomials of [2] are obtained by taking
% linear combinations of the Foster polynomials with
% random coefficients from the standard normal distribution:
ranpoly = @(n) foster(2:n)*randn(n-1,1);

%%
% Here, for example, are random polynomials of degrees 20, 100, and 500:
for n = 20*5.^(0:2)
    rng(2)
    plot(ranpoly(n)), ylim([-1 1]), grid on, snapnow
end

%%
% Note that, as usual for polynomials, these shapes are in no sense
% translation-invariant, having faster oscillations near the endpoints
% than in the middle.  Still, it is proved in [2] that they approach
% Brownian bridge as $n\to \infty$.  Many other interesting properties
% are also developed in [2], concerning polynomial moments for example,
% as well as applications to the numerical solution of stochastic
% differential equations.

%%
% How do random polynomials of a finite degree $n$ differ
% from true Brownian bridge?  
% It has been shown in [3] that the variance of the difference
% converges to zero at the rate $O(1/n)$
% and approaches a semicircle in profile.  An explicit formula
% involves $t-t^2$ plus the sum of the squares of the Foster
% polynomials.  Thus for example, here is the distribution
% for degree $n=20$:
s = chebfun('t-t^2',[0,1]);
for k = 2:20
  s = s - foster(k).^2;
end
plot(s), grid on

%%
% Beautiful!
%% References
% [1] S. Filip, A. Javeed, and L. N. Trefethen,
% Smooth random runctions, random ODEs, and Gaussian processes,
% _SIAM Review_ 61 (2019), 185--205.
%
% [2] J. Foster, T. Lyons, and H. Oberhauser,
% An optimal polynomial approximation of Brownian motion,
% arXiv:1904.06998, 2019.
%
% [3] K. Habermann, A semicircle law and decorrelation
% phenomena for iterated Kolmogorov loops,
% arXiv:1904.11484, 2019.
%
% [4] J.-P. Kahane, _Some Random Series of Functions_, 
% 2nd ed., Cambridge, 1985.
