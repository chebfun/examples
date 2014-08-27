%% The tiger's tail
% Nick Trefethen, August 2014

%%
% (Chebfun example roots/Tiger.m)
% [Tags: #ROOTS]

%%
% My essay "Six myths of polynomial interpolation and
% quadrature", reproduced as an Appendix in [1], closes
% with an example that reminds one of a tiger's tail.
% Here with a few modifications is that example:
x = chebfun('x',[-2 1]);
LW = 'linewidth'; MS = 'markersize';
CO = 'color'; orange = [1 .5 .25];
f = exp(.5*x).*(sin(5*x) + sin(101*x));
roundf = round(f);
r = roots(f-roundf);
hold off, plot(f,LW,2,CO,orange), hold on
ylim([-4 3])
plot(r,f(r),'.k',MS,12), hold off

%%
% Let's look at what's going on here.  First of all a
% chebfun $f$ is constructed:
plot(f,LW,1.6,CO,orange)
ylim([-4 3])

%%
% Then another chebfun is constructed consisting of $f$ rounded
% to integers:
plot(roundf,LW,1.1,'-k','jumpline','k')
ylim([-4 3])

%%
% Superimposing the two curves yields a lot of intersections,
% which are computed by `roots`:
number_of_roots = length(r)
plot(f,LW,1,CO,orange), hold on
plot(roundf,LW,.9,'-k','jumpline','k')
plot(r,f(r),'.k',MS,8), hold off

%%
% Notice that a dot appears not only where $f$ is equal to
% an integer, but also where it is equal to a half-integer.
% Do you see why?

%%
% Reference:
%
% 1. L. N. Trefethen, _Approximation Theory and Approximation
% Practice_, SIAM, 2013.
