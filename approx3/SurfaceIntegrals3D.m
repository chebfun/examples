%% Integration over 2D surfaces in 3D
% Olivier S&egrave;te, June 2016

%%
% (Chebfun Example approx3/SurfaceIntegrals3D.m)

%%
% In this example we illustrate the integration of scalar functions and vector
% fields over 2D surfaces in 3D space with the |integral2| command.

%% 
% Assume we have a surface given by a parametrization
% $S : D = [a,b] \times [c,d] \to S(D)$,
% $(u,v) \mapsto S(u,v)$, represented by a 3-component chebfun2v (i.e.,
% each point in the 2D domain $D$ is mapped to a point in 3D).

%%
% Let $f = f(x,y,z)$ be a scalar function defined over a box that contains the
% surface (i.e., $S(D)$).  We can then integrate $f$ over the surface with the
% following surface integral:
% $$\int_S f \, dS = \int_D f(S(u,v)) \left\Vert \frac{\partial S}{\partial 
% u}(u,v) \times \frac{\partial S}{\partial v}(u,v) \right\Vert \, dudv. $$
% When $f$ is represented by a chebfun3 object, this integral can be computed
% with |integral2|.

%%
% For a vector field $F(x,y,z) = [F_1(x,y,z); F_2(x,y,z); F_3(x,y,z)]$ defined
% on a box around containing the surface, we can compute its flux through the
% surface by the flux integral
% $$\int_S F \cdot \vec{dS} = \int_D F(S(u,v)) \cdot \left( \frac{\partial
% S}{\partial u}(u,v) \times \frac{\partial S}{\partial v}(u,v) \right) \,
% dudv. $$
% This can be done with |integral2|.

%%
% Let us consider the scalar function $f$ and the vector field $F$:

format long
cheb.xyz
f = x.^2 + y.*z;
F = [x+y; x.*z + y; cos(z)];

%%
% As our first example, let us consider the rippled disk parametrized by

S = chebfun2v(@(r, t) r .* cos(t), @(r, t) r.*sin(t), @(r, t) cos(5*r), ...
    [0, 5, 0, 2*pi]);
surf(S), axis equal

%%
% To compute the flux of $F$ through $S$ we simply type
integral2(F, S)


%%
% As a second example, we take $S$ to be the lower half of the unit sphere,
% parametrized by

S = chebfun2v(@(phi, theta) sin(theta) .* cos(phi), ...
    @(phi, theta) sin(theta) .* sin(phi), @(phi, theta) cos(theta), ...
    [0, 2*pi, pi/2, pi]);
surf(S), axis equal

%%
% The integral of $f$ over this surface is
integral2(f, S)

%%
% and the flux of $F$ through the bowl is
integral2(F, S)

%%
% Of course we can also integrate over curves in Chebfun3.  This is done with
% |integral|.  See 
% http://www.chebfun.org/examples/fun/LineIntegral3D.html
% for integration of a scalar function over a curve and
% http://www.chebfun.org/examples/fun/GaussGreenStokes.html
% for integrals of vector fields along a curve.
