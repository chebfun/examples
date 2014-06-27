%% Fourier-based chebfuns
% Grady Wright, June 2014

%%
% [Tags: #periodic, #interpolation, #fourier]

%%
% One of the new features of Chebfun version 5 is the ability to create
% chebfuns of smooth periodic functions using Fourier series. This example
% introduces and demonstrates some of the functionality of this new tool.
LW = 'linewidth'; lw = 1.6; MS = 'MarkerSize'; ms = 10;

%% Construction and comparison
% Fourier-based chebfuns, or "fourfuns" as we like to refer to them, can be
% created with the use of the `'periodic'` flag in the chebfun constructor.
% For, example, the function $f(x) = \cos(8\sin(x))$ for $-\pi \leq x \leq
% \pi$ can be constructed as follows:
% 
dom = [-pi,pi];
f = chebfun(@(x) cos(8*sin(x)),dom,'periodic')
plot(f,LW,lw);

%%
% Here $f$ is represented to machine precision using a Fourier interpolant
% rather than a Chebyshev interpolant. The displayed information for $f$ above
% shows that it is of length 61, meaning that $f$ is resolved to machine
% precision using 61 samples, or $(61-1)/2=30$ (complex) Fourier modes. These
% coefficients can be displayed graphically by
plotcoeffs(f)

%%
% Since $f$ is smooth and periodic, a Fourier representation requires fewer
% terms than a Chebyshev representation of $f$ to reach machine precision.
% We can check this by constructing $f$ without the `'periodic'` flag:
f_cheby = chebfun(@(x) cos(8*sin(x)),dom)

%%
% The ratio of length of the Chebyshev series to the Fourier series should
% be approximately $\pi/2$ since the former has a resolution power of 
% $\pi$ points per wavelength and the latter of 2 points per wavelength. 
% We can check this numerically as
ratio = length(f_cheby)/length(f)
theoretical = pi/2

%%
% Trying to construct a fourfun from a non-periodic or non-smooth function
% will typically result in a warning being issued and an "unhappy" fourfun,
% as illustrated for the unit step function below:
f = chebfun(@(x) 0.5*(1+sign(x)),dom,'periodic')
plot(f,LW,lw);

%%
% The length of $f$ is 65536, which is the maximum number of samples used
% in the construction process to try to resolve $f$. The famous Gibbs'
% phenomenon can be seen near the discontinuity in the plot of $f$. Chebfun
% can be used to represent this function in non-periodic mode (i.e. using
% Chebyshev series) with the option of `splitting on`:
f = chebfun(@(x) 0.5*(1+sign(x)),dom,'splitting','on')

%%
% Splitting is not an option for fourfuns.

%% Basic operations
% Many Chebfun operations can also be applied directly to a fourfun. 
% Some of these basic operations are illustrated in the examples below.

%%
% Addition, subtraction, multiplication, division, and function composition
% can all be directly applied to a fourfun.  However one should be aware that
% operation should result in a smooth and periodic function. The following
% example illustrates some of these operations:
g = chebfun(@(x) sin(x),dom,'periodic');
f = tanh(cos(1+2*g).^2)-0.5
plot(f, LW, lw)

%%
% The max, min, and roots of $f$ can be computed by
[maxf,xmaxf] = max(f);
[minf,xminf] = min(f);
rootsf = roots(f);
maxf
minf
rootsf

%%
% These can be visualized as
plot(f, LW, lw), hold on
plot(xmaxf,maxf,'gs',xminf,minf,'md',rootsf,0*rootsf,'ro',MS,ms)
legend('f','max f','min f','zeros f','location','southwest')
hold off;

%%
% The derivative of $f$ is computed using `diff`:
df = diff(f);
plot(df, LW, lw)

%%
% and the definite integral is computed using `sum`:
intf = sum(f)

%%
% Complex-valued fourfuns are also possible. For example:
f = chebfun(@(x) 1i*(13*cos(x)-5*cos(2*x)-2*cos(3*x)-cos(4*x)) + ...
                 16*sin(x).^3, dom, 'periodic')
plot(f, LW, lw), axis equal

%%
% The area enclosed by this curve can be computed as
area_heart = abs(sum(real(f).*diff(imag(f))))

%%
% According to [1], the true area enclosed is $180\pi$. The relative error
% in the computation above is then
err = (area_heart - 180*pi)/(180*pi)

%%
% The convolution of two smooth periodic functions can be computed using
% the `circconv` (circular convolution) function. The example below 
% demonstrates this function in combination with the additional feature 
% that allows fourfuns to be constructed from function values. The latter
% is demonstrated first:
rng('default'), rng(0);
n = 201;
x = fourpts(n);
func_vals = exp(sin(2*pi*x)) + 0.05*randn(n,1);
f = chebfun(func_vals,dom,'periodic')

%%
% Here $f$ interpolates the noisy `func_vals` at 201 equally spaced points
% from $[-\pi,\pi)$ using the Fourier basis. The high frequencies in this
% function can be smoothed by convolving it with a mollifier, in this case
% a (normalized) Gaussian with variance 0.1.
sigma = 0.1;
g = chebfun(@(x) 1/(sigma*sqrt(2*pi))*exp(-0.5*(x/sigma).^2),dom,'periodic');

%%
% Note that the resulting respresentation of $g$ is actually the periodic 
% extension of the Gaussian over $[-\pi,\pi]$.  The convolution of $f$ and
% $g$ is computed and visualized using
h = circconv(f,g);
plot(g,'b',LW,lw), hold on
plot(f,'r',LW,lw), plot(h,'k',LW,lw)
legend('Molifier g','Noisy function f','Smoothed function h');
hold off;

%% References
%
% 1. Mathworld Heart Curve: http://mathworld.wolfram.com/HeartCurve.html
