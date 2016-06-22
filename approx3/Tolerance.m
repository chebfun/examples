%% Loosening the Chebfun3 tolerance
% Nick Trefethen, June 2016

%%
% (Chebfun3 Example: Fun/HelloWorld3.m)

%% 1. Tolerances in Chebfun
% Chebfun's default tolerance is machine precision in 1D, 2D, and
% 3D:
chebfuneps
chebfun2eps
chebfun3eps

%%
% In 1D, there is usually
% not much to be gained by loosening the tolerance (unless you are
% working with noisy functions), and we have long recommended that users
% leave `chebfuneps` at its factory value.  (There is an
% FAQ question at `www.chebfun.org` on this topic.)  In 2D and especially 3D,
% however, loosening the tolerance is often worthwhile.
% This is discussed in Section 18.10 of the _Chebfun Guide_.

%%
% The reason the default tolerance is machine precision is
% that accurate results are often easily achievable.  For example, suppose
% we want to compute the triple integral
% $$ I = \int_{-1}^1 \int_{-1}^1 \int_{-1}^1 \exp(\sin(xyz + \exp(xyz))) dz dy dx . $$
% We could do it like this,
tic
f = chebfun3(@(x,y,z) exp(sin(x.*y.*z + exp(x.*y.*z))));
format long
I = sum3(f)
toc

%%
% We could also do it like this:
tic
cheb.xyz
f = exp(sin(x.*y.*z + exp(x.*y.*z)));
I = sum3(f)
toc

%%
% These results are quite satisfactory,
% because this chebfun3 is of only medium complexity:
f
[m,n,p] = length(f)

%%
%% 2. Slowdown for complicated functions
% On the other hand, if we make the function more complicated,
% things slow down:
tic
g = chebfun3(@(x,y,z) exp(sin(10*x.*y.*z + exp(x.*y.*z))));
I = sum3(g)
toc

%%
% Here are the parameters of the more complicated function:
g
[m,n,p] = length(g)

%% 3. Speedup if the tolerance is loosened
% A considerable speedup can often be achieved by working
% with a looser tolerance.  One way to construct the chebfun3
% with tolerance $10^{-8}$ is like this:
tic
g = chebfun3(@(x,y,z) exp(sin(10*x.*y.*z + exp(x.*y.*z))),'eps',1e-8);
I = sum3(g)
toc

%%
% Note that the value of $I$ agrees with the previous result
% to quite a few digits.  Here is the newly constructed chebfun3:
g
[m,n,p] = length(g)

%%
% Sometimes one wants to loosen the tolerance globally, e.g. if there
% will be further computations, like this:
chebfun3eps 1e-8
tic
cheb.xyz
g = exp(sin(10*x.*y.*z + exp(x.*y.*z)));
I = sum3(g)
toc

%%
% Let's try just four digits:
chebfun3eps 1e-4
tic
g = exp(sin(10*x.*y.*z + exp(x.*y.*z)));
I = sum3(g)
toc

%% 
% The computed integral still has several correct digits.

%%
% As good citizens, we now return the tolerance to its
% factory value:
chebfun3eps factory
