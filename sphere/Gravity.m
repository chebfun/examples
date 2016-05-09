%% Gravitational attraction to a sphere
% Nick Trefethen, May 2016

%%
% (Chebfun example sphere/Gravity.m)
% [Tags: #spherefun]

%% 1. A little geometry
% Here's a vector $X$ whose 2-norm is 1.5:
X = [-1  -1.1, -.2];
norm(X)

%%
% Let's use the vector-valued part of
% Spherefun to define the field of vector distances
% between points on the unit sphere and $X$, which we can
% think of as a spacecraft in orbit:
d = spherefunv(@(x,y,z) X(1)-x, @(x,y,z) X(2)-y, @(x,y,z) X(3)-z);

%%
% Here is the scalar function representing the _distance_ $r$ between
% $X$ and points on the sphere:
r = sqrt(d(1).^2 + d(2).^2 + d(3).^2); 

%%
% The minimum distance of $X$ to the sphere is of course $0.5%:
min_distance = min2(r)

%%
% Similarly, the maximum is $2.5$
max_distance = max2(r)

%% 
% Here is a contour plot of $r$ on the sphere,
% together with a red dot showing our little spacecraft.
contour(r,.6:.1:2,'k')
hold on, plot3(X(1),X(2),X(3),'.r','markersize',25), hold off
view(-10,35)

%% 2. Inverse-square gravitational force
% A great discovery of Newton (or was it Hooke?) is that the gravitational
% forces associated with a sphere of uniform mass distribution are
% the same as if all the mass were concentrated at the center.
% Accordingly, we know that if a unit mass is spread around the
% sphere and the spacecraft also has unit mass, then the 
% inverse-square attraction between them should be $(1.5)^{-2}$,
% since the distance from the spacecraft to the center of the sphere is $1.5%:
force_exact = 1/1.5^2

%%
% Let's confirm Newton's observation by computing the net gravitational
% force as an integral over the surface of the sphere.
% Since the area of the sphere is $4\pi$, the density of a uniformly
% distributed mass is
rho = 1/(4*pi)

%%
% That gives us the following component of the force at each point,
% in the direction of $X$:
Xnormalized = X/norm(X);
force_function = rho*(Xnormalized*d)./r.^3;

%%
% Summing, we get the right answer:
force = sum2(force_function)
