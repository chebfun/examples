%% The Gibbs phenomenon in 2D
% Andre Uschmajew and Nick Trefethen, February 2017

%%
% (Chebfun example approx2/Gibbs2D.m)

%% 1. Chebyshev 2D Gibbs effect
% Here is an illustration of the Gibbs phenomenon in 2D:
A = zeros(100); A(40:61,40:61) = 1;
p = chebfun2(A); plot(p), view(-110,60)

%%
% What's going on is that Chebfun has constructed a bivariate
% polynomial interpolant $p(x,y)$ to data at $100^2$
% 2D Chebshev points, zero on most of the domain and 1 on a
% square in the middle.

%%
% How big is the oscillation?
min2(p), max2(p)

%%
% This is big! -- about twice what we are used to with a 
% 1D Gibbs effect:
a = zeros(100,1); a(40:61) = 1; 
p1 = chebfun(a); max(p1)

%%
% The reason is that the overshoot at the corner is especially
% large, as we can see by zooming in:
pzoom = p{0,.5,0,.5}; plot(pzoom), view(-100,20)

%% 2. Fourier 2D Gibbs effect
% A Fourier analogue can be produced by including the
% 'periodic' flag:
t = chebfun2(A,'periodic'); plot(t), view(-110,60)

%%
% The extrema are similar:
min2(t), max2(t)

%%
% For a summary of numbers associated with various versions
% of the Gibbs phenomenon, see Chapter 9 of [1].

%% Reference
%
% 1. L. N. Trefethen, _Approximation Theory and Approximation Practice_,
% SIAM, 2013.
