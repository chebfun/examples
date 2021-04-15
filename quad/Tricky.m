%% Some tricky integrals
% Fredrik Johansson and Nick Trefethen, March 2018

%%
% (Chebfun example quad/Tricky.m)
% [Tags: #quadrature]

%% 0. Introduction
% FJ gave a talk at ENS Lyon today with a number of examples
% in it that intrigued LNT.  Here we play with some of those
% examples in Chebfun.  Now Chebfun is just numerical, with
% no guarantees of accuracy, which means it is not a competitor
% for FJ's Arb method of rigorous quadrature [1].  The point of this
% example is only to see how Chebfun does, unrigorously in
% floating point arithmetic,
% on some challenging examples people have cooked up over the years.

%%
% We find that in all but one of these examples, Chebfun very nicely
% gets a highly accurate answer.  Example 3, however, shows
% difficulties with a function with a lot of discontinuities.

%% 1. Three spikes
% This problem comes from R. Cranley and T. N. L. Patterson, On
% the automatic numerical evaluation of definite integrals,
% The Computer Journal 14 (1971), 189-198.  It also appeared in
% D. K. Kahaner, Comparison of numerical quadrature formulas, in
% J. R. Rice, ed., Mathematical Software, Academic Press, 1971,
% 229-259.
% See also the earlier Chebfun example
% www.chebfun.org/examples/quad/SpikeIntegral.html.
format long
ff = @(x) 1/cosh(10*(x-.2))^2 + 1/cosh(100*(x-.4))^4 + 1/cosh(1000*(x-.6))^6;
Iexact = 0.210802735500549277
tic, f = chebfun(ff,[0 1]), I = sum(f); toc
plot(f)

%%
% Note that turning on splitting doesn't make much difference to speed.
tic, f = chebfun(ff,[0 1],'splitting','on'), I = sum(f); toc

%% 2. Violent oscillation
% This problem comes from S. M. Rump, Verification methods:
% rigorous results using floating-point arithmetic,
% Acta Numerica, 19 (2010), 287-449.  It is also discussed
% in W. Tucker, Validated Numerics: A Short Introduction to
% Rigorous Computations, Princeton U. Press, 2011.

ff = @(x) sin(x+exp(x));
Iexact = 0.34740017265724780787
tic, f = chebfun(ff,[0 8]), I = sum(f), toc
plot(f)

%% 3. Violent oscillation with 2979 discontinuities
% If we try this with default parameters, we get a warning related
% to Chebfun's default preference values
% splitMaxLength = 6000.
ff = @(x) (exp(x)-floor(exp(x)))*sin(x+exp(x));
Iexact = 0.098651704478365206119
tic, f = chebfun(ff,[0 8],'splitting','on'); I = sum(f), toc

%%
% By increasing splitMaxLength greatly, we can get an answer
% but it's outrageously slow *and* it's only accurate to 6 digits.
% Specifically, the commands
%
%   I = sum(chebfun(ff,[0 8],'splitting','on','splitMaxLength',1e6))
%
% give the result I = 0.0986522613...  after 72 seconds.  We are
% well
% aware that Chebfun is slow for problems with many discontinuities,
% but for it to lose ten digits of accuracy looks like a bug
% somewhere.  (I have confirmed by adding up the pieces by hand that
% the ``exact'' answer is correct.)

%%
% It would be good to investigate why this integral is giving such
% trouble.

%% 4. Error function
% This problem comes from Silviu Filip.
ff = @(x) exp(-x)*erf(sqrt(1250)*x+1.5);
Iexact = NaN
tic, f = chebfun(ff), I = sum(f), toc
plot(f)

%% 5. Airy function
% This problem comes from FJ.
ff = @(x) exp(-x)*airy(-x);
Iexact = 0.378751605379086535
tic, f = chebfun(ff,[0 inf]), I = sum(f), toc
plot(f)

%%
% We compare this with the result on a sufficiently large
% finite interval:
tic, f = chebfun(ff,[0 40]), I = sum(f), toc

%% 6. Absolute value of polynomial
% This problem was posed by Harald Helfgott on MathOverflow,
% https://mathoverflow.net/questions/123677/rigorous-numerical-integration.
% See also A. Mahboubi, G. Melquiond, and T. Sibut-Pinote,
% Formally verified approximations of definite integrals,
% International Conference on Interactive Theorem Proving,
% Spring, 2016, 274-289.
ff = @(x) abs(x^4+10*x^3+19*x^2-6*x-6)*exp(x);
Iexact = 11.1473105500571397339
tic, f = chebfun(ff,[0 1],'splitting','on'), I = sum(f), toc
plot(f)

%% 7. A ceiling function
% This is an integral representation of the triangular
% sum that Gauss famously figured out as a schoolboy.
% See B. Hayes, Gauss's day of reckoning, American Scientist
% 94 (2006), 200-205.
ff = @(x) ceil(x);
Iexact = 5050
tic, f = chebfun(ff,[0 100],'splitting','on'); I = sum(f), toc
plot(f)

%% 8. Another non-smooth function
ff = @(x) (x-floor(x)-1/2)*max(sin(x),cos(x));
Iexact = -0.14281864202632808376
tic, f = chebfun(ff,[0 10],'splitting','on'); I = sum(f), toc
plot(f)

%% 9. From Brisebarre and Joldes
% This example comes from Mioara Joldes,
% Rigorous Polynomial Approximations and Applications,
% PhD thesis, ENS Lyon, 2011.  The source of the integral is
% C.-Y. Chen, Computing interval enclosures for definite
% integrals by application of triple adaptive strategies,
% Computing, 78 (2006), 81-99.
% It requires a Chebfun of length more than a million.
ff = @(x) sin((0.001+(1-x)^2)^(-3/2));
Iexact = 0.74997436852719477011
tic, f = chebfun(ff,[0 3],'maxLength',1e7), I = sum(f), toc
plot(f)

%%
% It is interesting to note the near-zero region near x=1 here.  That's
% incorrect, but has negligible effect on the integral.  

%%
% Splitting on is hard work too, but at least it eventually
% gets the right answer, and this time with a more convincing plot.
tic, f = chebfun(ff,[0 3],'splitting','on','splitMaxLength',1e6);
I = sum(f), toc
plot(f)

%% Reference
% Not much has been done on rigorous extended precision arithmetic,
% but FJ's Arb library is a contribution in this area:
% 
% [1] F. Johansson, Numerical integration in
% arbitrary-precision ball arithmetic, arXiv:1802.07942 and
% _International Congress on Mathematical Software_, Springer, Cham, 2018.
