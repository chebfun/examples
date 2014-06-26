%% Integrals over closed contours using periodic chebfuns
% Mohsin Javed, June 2014

%%
% (Chebfun example complex/ClosedContours.m)
% [Tags: #contour, #contourintegral, #periodic]

%%
%
LW = 'linewidth'; lw = 1.5;

%%
% In this example, we compute a few integrals over closed contours in the complex
% plane using periodic chebfuns.
%%
% Consider a smooth and closed contour $\Gamma$ in the complex plane and let
% us say we want to compute
% $$ \int_{\Gamma} f(z) dz. $$
% If we parametrize $\Gamma$ using a real varaible, say $t$ , then since the
% contour is closed, the intgrand becomes periodic in $t$ and we get
% $$ \int_{\Gamma} f(z) dz = \int_{a}^{b} f(z(t)) z'(t) dt. $$
%%
% All this can be done very efficiently in Chebfun. Thanks to the underlying
% Fourier-technology which has been integrated with Chebfun's usual Chebyshev
% technology.

%%
% Here is a simple example. Consider the function:
f = @(z) (1-2*z)./(z.*(z-1).*(z-3));

%%
% Say we want to integrate this funciton on a circle of radius $2$. To do this
% in Chebfun's periodic mode, we first parametrize the circle:
z = chebfun(@(t) 2*exp(2*pi*1i*t), [0, 1], 'periodic');
plot(z, LW, lw)
ylim( [-2.5, 2.5] ), axis equal
title('Contour of Integration')
%%
% The integrand is then constructed by a simple composition:
I = f(z)
%%
% This is how the real and imaginary part of the integrand looks on the contour:
subplot(1, 2, 1)
plot(real(I), LW, lw)
title('real part of the function')
subplot(1, 2, 2)
plot(imag(I), LW, lw)
title('Imaginary part of the function')
%%
% To compute the integral, we recall that
% $$ \int_{|z|=2} f(z) dz = \int_{0}^{1} f(z(t)) z'(t) dt. $$
% We theefore, first compute $z'(t)$:
dz = diff(z);
%%
% Computing the integral now can not be easier:
s = sum(I.*dz)

%%
% The true answer is $5 \pi i/3$, and we see that Chebfun has done a very good
% job:
norm(s - 5/3*pi*1i)

%%
% Here is another example. Consider the well known function
f = @(z) sin(z)./z;
%%
% This analytic function has a remobable singularity at the origin. Therefore,
% the integral of the function on any closed contour should be zero according to
% Cauchy's theorem.
z = chebfun(@(t) exp(2*pi*1i*t), [0, 1], 'periodic');
I = f(z);
dz = diff(z);
%%
% Again, we plot the real and imaginary part of the integrand on the unit
% circle:
subplot(1, 2, 1)
plot(real(I), LW, lw)
title('real part of the function')
subplot(1, 2, 2)
plot(imag(I), LW, lw)
title('Imaginary part of the function')

%%
% And here is the integral, which is numerically zero:
s = sum(I.*dz)

%%
% As our final example, we pick a function with an essential singularity at the
% origin and compute its integral on the unit circle.
f = @(z) exp(1./z).*sin(1./z);
z = chebfun(@(t) exp(2*pi*1i*t), [0, 1], 'periodic');
I = f(z);
dz = diff(z);
s = sum(I.*dz)