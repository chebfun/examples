%% Simple computations with probability distributions
% Mark Richardson, May 2011

%%
% (Chebfun example stats/Expectations.m)

%%
% In this example, we use Chebfun to solve some probability distribution
% problems from [1].

%% 1. Expectation of a random variable
%
% We use Problem 3.4 from p. 86 of [1] to motivate this example.
%
% Suppose a continuous random variable $X$ has the probability density function
%
% $$ f(x) = 2e^{-2x},~~  x \ge 0,  \qquad
% f(x) = 0,  ~~    x < 0. $$
%
% What are: (a) $E(X)$ and (b) $E(X^2)$?

%%
% (a) In order to compute the expectation $E(X)$, we first define a
% chebfun for $X$.  This can be done over the semi-infinite interval
% $[0,\infty)$, but the resulting integrals lose a few digits of
% precision.  Instead, since $\exp(-x)$ decreases so quickly, we
% take the interval to be $[0,40]$. 
x = chebfun('x',[0 40]);

%%
% Next we approximate the density function.
f = 2*exp(-2*x);
figure('Position',[100 200 520 180])
LW = 'linewidth'; lw = 1.6;
plot(f,LW,lw), grid on
ylim([-0.2 2.2])
xlabel('x'), ylabel('f(x)','rotation',0)

%%
% If $f$ is a density function, its integral should be $1$, and we find
% that this is indeed the case to within rounding errors.
sum(f)

%%
% The expectation of a continuous random variable is defined as the
% integral over of $xf(x)$.
xf = x.*f;
plot(xf,LW,lw), grid on
ylim([-0.05 0.4])
xlabel('x'), ylabel(sprintf('x f(x)\n'),'rotation',0)

%%
% We can use the chebfun command `sum` to compute this integral
% The correct answer in this case is $1/2$.
format long
sum(xf)

%%
% b) For $E(X^2)$, the answer is again $1/2$ and we compute
% this in the same way as before.
xxf = x.^2.*f;
plot(xxf,LW,lw), grid on
ylim([-0.03 0.31])
xlabel('x'), ylabel('x^2 f(x)','rotation',0)
%%
sum(xxf)

%% 2. Mean, median and mode of a probability distribution
% This example is motivated by problem 3.33 on p. 94 of [1].

%%
% The probability density function of a continuous random variable $X$ is
%
% $$ g(x) = 4x(9-x^2)/81, ~~ 0\le x\le 3, $$
%
% and zero otherwise. Find: a) the mean,
% b) the median, and c) the mode.

%%
% First, we define an appropriate Chebfun variable and the pdf:
x = chebfun('x',[0 3]);
g = 4*x.*(9-x.^2)/81;
plot(g,LW,lw), grid on
ylim([-0.01 0.61])
xlabel('x'), ylabel('g(x)','rotation',0)

%%
% a) Computing the mean is simply a matter of computing the expectation as
% in the previous question. The exact answer is $1.6$ and this is what we
% find using Chebfun.
mean = sum(x.*g)

%%
% b) The median is the value $a$ for which $P(X\le a) = 1/2$.
% In order to solve
% this problem we need to work with the cumulative distribution function,
% which is simply the indefinite integral of the probability density. This
% can be computed with the chebfun command `cumsum`.
G = cumsum(g);
plot(G,LW,lw), grid on
xlabel('x'), ylabel(sprintf('G(x)\n'),'rotation',0)

%%
% Note again that as we would expect for any pdf, the integral is $1$. 
% Here is the median $a$:
median = roots(G-0.5)
median_exact = sqrt(9-9*sqrt(2)/2)

%%
% c) For the mode, we are looking for the position of the global
% maximum of the probability distribution. This is easily computed
% with the Chebfun command `max`.
[gmax,mode] = max(g);
display(mode)

%%
% Again, this matches the exact result
mode_exact = sqrt(3)

%%
% Here is a graph showing the three computed values:
plot(g,LW,lw), grid on, hold on
plot([mean mean],[0 g(mean)],'-r',LW,lw)
plot([median median],[0 g(median)],'-m',LW,lw)
plot([mode mode],[0 g(mode)],'-k',LW,lw)
text(0.2,0.55,sprintf('mean   = %1.2f',mean),'color','r')
text(1.2,0.55,sprintf('median = %1.2f',median),'color','m')
text(2.2,0.55,sprintf('mode   = %1.2f',mode),'color','k')
hold off, ylim([-0.01 0.61])
xlabel('x'), ylabel('g(x)','rotation',0)

%% Reference
%
% 1. M. Spiegel, J. Schiller, and R. Srinivasan, _Schaum's Outlines:
%    Probability and Statistics_, 3rd. ed., 2009.
