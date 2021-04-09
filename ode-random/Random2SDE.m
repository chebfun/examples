%% From random functions to SDEs
% Nick Trefethen and Abdul-Lateef Haji-Ali, May 2017

%%
% (Chebfun example ode-random/Random2SDE.m)

%%
% The new Chebfun release v5.7.0 brings with it a new
% set of capabilities: smooth random functions.
% The commands are |randnfun|, |randnfun2|, |randnfunsphere|, and 
% |randnfundisk|.  (There is not yet
% a |randnfun3|.)  Each of these returns a band-limited
% function defined by a Fourier series with independent random
% coefficients; a parameter $\lambda$ specifies the
% minimal wavelength.  For information about the mathematics, see [1] and
% also chapter 12 of [2].
help randnfun

%%
% These commands make it easy to compute sample paths of random
% ODEs in Chebfun, that is, of ODEs defined by smooth random coefficients.
% As $\lambda \to 0$, we approach the limit of stochastic DEs (SDEs)
% in their Stratonovich (as opposed to It&ocirc;) formulation.
% Thus Chebfun provides an easy window into
% the phenomena that make SDEs so fascinating and so important.  Note that
% what Chebfun solves is a random ODE
% (based on band-limited randomness), not an SDE (based on band-unlimited
% randomness, i.e., white noise, a notion made precise through the 
% formulation of a Wiener process, also known as Brownian motion).
% In the limit $\lambda\to 0$ they are the same (though Chebfun becomes
% quite an inefficient tool if you try to get too close to that limit).

%%
% For ODE-related studies, one should always call |randnfun| and its
% cousins with the flag |'big'|.  This multiplies the random function
% by $(\lambda/2)^{-1/2}$, meaning that its amplitude grows without bound
% as $\lambda\to 0$.  This is what is needed for random ODEs to approximate SDEs.

%%
% Here we give just the simplest example.  If $f$ is a normalized random
% function, then the equation $$ u' = f $$
% has the indefinite integral of $f$ as its solution.  We call this a
% "smooth random walk".  If $\lambda$ is small enough, it looks to the eye
% like true Brownian motion, and as $\lambda \to 0$, that is what it approaches.
% For finite $\lambda$ there are no mathematical technicalities to worry about;
% it is simply an ODE.  Precise statements about the limit $\lambda\to 0$
% require care, however, and stochastic analysts would write the equation
% above in a very different form,
% $$ dX_t = d W_t. $$

%%
% Here are three smooth random walk sample paths for $t\in [0,1]$ with $u(0) = 0$
% and $\lambda = 0.001$.
tic
rng(0)
u = randnfun(0.001,[0 1],3,'big');
plot(cumsum(u))
grid on, ylim([-2 2])
xlabel t, ylabel u

%%
total_time_in_seconds = toc

%%
% References:
%
% [1] S. Filip, A. Javeed, and L. N. Trefethen,
% Smooth random functions, random ODEs, and Gaussian processes,
% _SIAM Review_, 61 (2019), 185--205.
%
% [2] L. N. Trefethen, _Approximation Theory and Approximation Practice,
% Extended Edition_, SIAM, 2019.
