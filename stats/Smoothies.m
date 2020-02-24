%% Smoothies: nowhere analytic functions
% Nick Trefethen, February 2020

%%
% (Chebfun example stats/Smoothies.m)
% [Tags: smoothie, #randnfun]

%%
% Since Weierstrass in the 19th century, we have known that there are
% functions that are continuous and yet nowhere differentiable.
% Such functions are associated with Fourier or Chebyshev series that
% are lacunary or random.  A beautiful example is a Brownian path,
% which can be defined by a Fourier series with random coefficients
% of magnitudes decreasing inverse-linearly [1,2].

%%
% At the other end of the smoothness spectrum, in a common room
% discussion in Oxford we found ourselves wondering, what about
% functions that are 
% $C^\infty$, i.e. infinitely differentiable, but nowhere analytic?
% It is well known how to cook up a $C^\infty$ function that is
% nonanalytic at a single point.  The standard example is
% $f(x) = \exp(-1/x^2)$ for $x\in [-1,1]$, which has zero derivatives of 
% all orders at $x=0$ and thus a Taylor series $0 + 0x + 0x^2 + \cdots,$ but
% with the Taylor series obviously not converging to $f$.  Can we
% upgrade this to a function which has Taylor series everywhere in
% $[-1,1]$, yet with none of them converging to the right limit?

%%
% Random Fourier series again give an elegant solution, and in 
% Chebfun, with apologies to the fruit beverage industry,
% this has been implemented by the new command |smoothie|.
% Here is a smoothie on the default interval $[-1,1]$:
rng(1)
f = smoothie;
plot(f)
ylim([-4 4])

%%
% To see the idea, we can look at the magnitudes of
% the Chebyshev coefficients:
plotcoeffs(f,'color','k')

%%
% These are random numbers of amplitudes decreasing
% faster than the reciprocal of any polynomial,
% but slower than exponentially.  This is enough to guarantee that
% (with probability 1) $f$ is $C^\infty$ but nowhere analytic.
% To be a little more precise, the construction of a smoothie closely
% follows that of a smooth random function, implemented in
% the command |randnfun|.  To ensure that the statistical properties
% are translation-invariant rather than changing as $x$ approaches
% $\pm 1$, the function is actually constructed via a random
% Fourier rather than Chebyshev series (hence periodic) on a longer interval
% than $[-1,1]$, which is then restricted to $[-1,1]$.  
% The Fourier coefficients decrease root-exponentially in amplitude,
% that is, at a rate $C^{-\sqrt n}$ with $C>1$.

%%
% To get a periodic smoothie, one can use the |'trig'| flag:
ftrig = smoothie('trig');
plot(ftrig), snapnow
plotcoeffs(ftrig,'color','k')

%%
% For all the details, take
% a look at the (quite simple) code |smoothie.m|.
% Another option is a complex smoothie,
fcomplex = smoothie('complex');
plot(fcomplex,'m'), ylim([-1.8 1.8]), axis equal

%%
% What do the Taylor series of these functions look like?  To give the
% idea, we plot the first and second derivatives of the original
% function $f$ at the beginning of this example:
subplot(2,1,1), plot(diff(f)), ylim([-80 80])
subplot(2,1,2), plot(diff(f,2)), ylim([-8000 8000])

%%
% These are smooth functions, but they are rapidly
% getting bigger in amplitude as well as having visible structure
% on smaller and smaller space scales.  At any point $x\in [-1,1]$, the Taylor
% series of $f$ will be well defined, but with coefficients 
% growing too fast for a positive radius of convergence.

%%
% References:
%
% [1] S. Filip, A. Javeed, and L. N. Trefethen,
% Smooth random functions, random ODEs, and Gaussian processes,
% _SIAM Review_ 61 (2019), 185--205.
%
% [2] J.-P. Kahane, _Some Random Series of Functions_, 
% 2nd ed., Cambridge, 1985.
