%% Integrating Tj(x)*Tk(y) over the unit disk
% Mikael Slevinsky, Nick Trefethen, and Klaus Wang, May 2016

%% 1. Numerical integration over the disk
% In studying cubature formulas, we needed
% to compute the integrals of products of Chebyshev polynomials
% $T_j(x) T_k(y)$ over the unit disk, like this one:
T8 = chebpoly(8); T16 = chebpoly(16);
s = linspace(-1,1,160); [xx,yy] = meshgrid(s,s);
ff = T8(xx).*T16(yy);
ff(xx.^2+yy.^2>1) = NaN;
contour(s,s,ff), axis equal off
hold on, plot(chebfun('exp(1i*x)',[0 2*pi]),'k'), hold off

%%
% If $j$ or $k$ is
% odd, the integral is zero, but even if they are both 
% even, to our surprise, we found that the
% integrals are still usually zero!  In fact, 
% a nonzero result only
% shows up if $j$ and $k$ differ by 0 or 2.
% Thus the function plotted above, for example,
% has integral exactly zero over the disk.  This is not obvious.

%%
% Speaking in general, suppose we want to integrate a
% smooth function $f(r,t)$ numerically over the unit disk, where $r$ is
% radius and $t$ is angle.  Soon Diskfun will be available for
% such problems, but here, we use standard Chebfun.  Let's
% assume that $f$ is defined by a function that accepts a scalar
% $r$ and a vector $t$, like this one:
f = @(r,t) 1 + 0*t;

%%
% Here is a numerical confirmation that the integral of $1$ is $\pi$:
format long
fr = @(r) r*sum(chebfun(@(t) f(r,t),[0,2*pi],'trig'));
I = sum(chebfun(@(r) fr(r),[0 1],'vectorize'))
Iexact = pi

%%
% For the function $f(r,t) = r^2$, the integral is $\pi/2$:
f = @(r,t) r^2 + 0*t;
fr = @(r) r*sum(chebfun(@(t) f(r,t),[0,2*pi],'trig'));
I = sum(chebfun(@(r) fr(r),[0 1],'vectorize'))
Iexact = pi/2

%%
% For the function $f(r,t) = r^2 \cos^2(t)$, the integral is $\pi/2$:
f = @(r,t) r^2*cos(t).^2;
fr = @(r) r*sum(chebfun(@(t) f(r,t),[0,2*pi],'trig'));
I = sum(chebfun(@(r) fr(r),[0 1],'vectorize'))
Iexact = pi/4

%% 2. Numerical integration of products of Chebyshev polynomials
% What about Chebyshev polynomials?  Here's a table.
tic
I = zeros(6);
format short
for j = 0:2:10
  Tj = chebpoly(j);
  for k = 0:2:j
    Tk = chebpoly(k); 
    f = @(r,t) Tk(r*cos(t)).*Tj(r*sin(t));
    fr = @(r) r*sum(chebfun(@(t) f(r,t),[0,2*pi],'trig'));
    I(1+j/2,1+k/2) = sum(chebfun(@(r) fr(r),[0 1],'vectorize'));
  end
end
I = I + tril(I,-1)'
time_elapsed_in_seconds = toc

%% 3. An analytic expression for the integrals

%% 4. Application to integration of general functions
% The results of the last section imply that it is very easy to
% compute the integral of a chebfun2 over a disk: just take the
% appropriate linear combination of its 
% bivariate Chebyshev coefficients.  For example [to be completed...
% and should we turn this into a new Chebfun2 command 
% called sumdisk?].

