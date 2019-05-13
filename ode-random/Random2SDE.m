%% From random functions to SDEs
% Nick Trefethen and Abdul-Lateef Haji-Ali, May 2017

%%
% (Chebfun example ode-random/Random2SDE.m)

%%
% The new Chebfun release v5.7.0 brings with it a new
% set of capabilities: smooth random functions, or as we
% will often simply call them, random functions.  The commands
% are |randnfun|, |randnfun2|, |randnfunsphere|, and 
% |randnfundisk|.  (As of this writing there is not yet
% a |randnfun3|.)  Each of these returns a band-limited
% function defined by a Fourier series with independent random
% coefficients; a parameter $\lambda$ specifies the
% minimal wavelength.  These commands have not yet been fully written up, but
% a few examples can be found in the Approximation (1D), Approximation (2D),
% and Statistics and Probability sections of these Examples, and
% more details can be found in the help texts and Matlab codes.
help randnfun

%%
% These commands make it easy to compute sample paths of random
% ODEs in Chebfun, that is, of ODEs defined by smooth random coefficients.
% As $\lambda \to 0$, we approach the limit of stochastic DEs (SDEs), typically
% in their Stratonovich (as opposed to It&ocirc;) formulation.  This too has
% not been written up.  Precise statements relating Chebfun's random
% functions to SDEs are not yet available, but we expect them to be
% developed in due course, building for example on the
% theory of Wong and Zakai [1].  But we don't need a fully developed theory
% to start exploring!  Chebfun provides an easy window into some of
% the phenomena that make SDEs so fascinating and so important.  Just
% remember that, properly speaking, what Chebfun solves is a random ODE
% (based on band-limited randomness), not a true SDE (based on band-unlimited
% randomness, i.e., white noise, a notion made precise through the 
% formulation of a Wiener process, also known as Brownian motion).

%%
% For studies of this kind, one should always call |randnfun| and its
% cousins with the flag |'big'|.  This multiplies the random function
% by $\lambda^{-1/2}$, meaning that its amplitude grows without bound
% as $\lambda\to 0$.  This is what is needed for random ODEs that
% are intended to approximate SDEs.

%%
% Here we give just the simplest example.  If $f$ is a normalized random
% function, then the equation
% $$ u' = f $$
% has the indefinite integral of $f$ as its solution.  We call this a
% "smooth random walk".  If $\lambda$ is small enough, it looks to the eye
% like true Brownian motion, and as $\lambda \to 0$, that is what it approaches.
% For finite $\lambda$ there are no mathematical technicalities to worry about;
% it is simply an ODE.  Precise statements about the limit $\lambda\to 0$
% require care, however, and stochastic analysts would write the equation
% above in a very different form,
% $$ dX_t = d W_t. $$

%%
% Here are three smooth random walk
% sample paths for $t\in [0,1]$ with $u(0) = 0$
% and $\lambda = 0.001$.
tic
rng(0)
u = randnfun(0.001,[0 1],3,'big');
plot(cumsum(u),'linewidth',2.5)
grid on, ylim([-1 1])
xlabel('t','fontsize',30), ylabel('u','fontsize',30)

%%
total_time_in_seconds = toc

%%
% Reference:
%
% [1] E. Wong and M. Zakai, On the convergence of ordinary
% integrals to stochastic integrals, _The Annals of Mathematical Statistics_
% 36.5 (1965), 1560-1564.
