%% Two electrons orbiting symmetrically about a nucleus
% Jeremy Fleury and Nick Trefethen, June 2016

%%
% (Chebfun example ode-nonlin/Electrons.m)

%% 1. Typical trajectories
% Here is a variation on the familiar $n$-body problem, suggested
% to us by Charlie Peskin of New York University.
% Suppose $n$ electrons of charge
% $-1$ are flying around a nucleus of infinite mass and
% charge $+n$.  What do the trajectories look like?  
% For $n=1$ it's trivial, just a circular orbit.  For $n\ge 2$ one
% sees all kinds of disordered and chaotic trajectories, and we
% hope to explore some of the
% possibilities in a future example.  In this example we consider
% just  a very special configuration with $n=2$, in which
% the two electrons are assumed to be exactly symmetrical about a line
% of reflection.  Here is a typical trajectory over a time interval
% of length 40.  
LW = 'linewidth'; MS = 'markersize'; lw = 1.2;
N = chebop(0,40);
N.op = @(t,z) diff(z,2) + 2*z./abs(z).^3 - 0.25i*imag(z)./imag(z).^3;
V = 1; N.lbc = [1i; V];
tic, z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off
set(gca,'xtick',-1:1,'ytick',-1:1)

%%
% For this computation we have used complex arithmetic for
% convenience, with the nucleus at the origin.  Because of the
% symmetry, only one particle needs to be tracked, so we
% have a scalar complex nonlinear second-order ODE initial-value
% problem:
% $$ z'' = {-2z\over |z|^3} + {z-\overline{z}\over |z-\overline{z}|^3}. $$
% For initial conditions in this example we take
% $$ z(0)=i, \quad z'(0) = V > 0. $$
% In the figure above, $V=1$.

%%
% Though it is not periodic, this orbit has a great deal of
% regularity, as we can see by plotting the $x$ component as
% a function of time.  
plot(x), xlabel t, ylabel('x(t)'), ylim([-1.5 1.5])

%%
% This is not a chaotic problem; it is more like quasiperiodic.
% As $t\to\infty$, the trajectory fills up
% a certain region in the $x\hbox{-}y$ plane.

%% 2. Energy
% The kinetic energy of this motion is $|z'|^2$, and the
% potential energy is $-4/|z| + 1/2\,\hbox{Im} z$.  Thus the total energy is
% $$ E = |z'|^2 -{4\over |z|} + {1\over 2\,\hbox{Im} z}, $$
% and this quantity is conserved.  For our initial value
% $V=1$, the energy is $E = -2.5$.  (We do not verify
% this by a Chebfun computation, which would be very slow because
% of the near-singularities when $\hbox{Im}z$ is near zero.)

%%
% Note that a particle at $z=\infty$ with zero velocity
% has energy $0$, and our initial condition will have energy $0$
% with this initial velocity:
% $$ V_{\hbox{crit}} = \sqrt{3.5} \approx 1.87. $$
% Sure enough, with further experiments one may verify that
% for $V>V_{\hbox{crit}}$, the
% trajectory flies off to infinity.  The reader may also find it
% interesting to explore "Pluto" trajectories just below this limit
% starting from values such as $V=1.85$ or $1.86$.

%% 3. Periodic trajectories
% For certain special initial velocities, the orbits are periodic.
% The simplest one, corresponding to $V\approx 1.446$, has
% the electrons simply swinging back and forth:
N.domain = [0 20];
V = 1.446; N.lbc = [1i; V];
z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off
set(gca,'xtick',-1:1,'ytick',-1:1)

%%
% A plot of $x$ values confirms the periodicity:
plot(x), xlabel t, ylabel x

%%
% Here is an estimate of the period $T$:
r = roots(x-.9*max(x)); r = r(deriv(x,r)>0); T = r(2)-r(1)

%%
% Here is another periodic solution:
V = 0.783; N.lbc = [1i; V];
z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off

%%
% Since the orbit is more complicated, it is not surprising that
% the period is longer:
set(gca,'xtick',-1:1,'ytick',-1:1)
r = roots(x-.9*max(x)); r = r(deriv(x,r)>0); T = r(2)-r(1)

%%
% And here a third periodic solution with a still longer period:
V = 1.17745; N.lbc = [1i; V];
z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off
set(gca,'xtick',-1:1,'ytick',-1:1)
r = roots(y-.9999999);
T = mean(r(end-1:end))

%%
% Readers trying these computations on their own may enjoy
% experimenting with the command |comet(z)|, which shows the planets
% flying about.

%% 4. Periodic representation
% A periodic function should have a trigonometric representation!
% We can compute such a representation with Chebfun's "trig" option.

%%
% First, let's make our estimate of the period a little more accurate.
% With the current estimate, instead of the expected value $1i$,
% we find this result after one cycle:
z(T)

%%
% An improved estimate of T is accordingly
T = T - real(z(T))/V

%%
% which gives a few more digits:
z(T)

%%
% Here we construct a trigonometric representation:
% using a loosened tolerance since the computed orbit and
% period are still far from machine precision.
% The match with the nonperiodic representation is good:
zT = z{0,T};
zTtrig = chebfun(zT,'trig','eps',1e-8);
norm(z{0,T}-zTtrig)

%%
% These are the absolute values of the Fourier coefficients:
plotcoeffs(zTtrig,'.k')

%%
% A plot of the velocities over the course of one orbit
% is interesting too:
plot(diff(zTtrig),'m')
title('Velocities z''(t)')

%% 5. Computing time
% Chebfun is wonderfully convenient, but it is not especially
% fast, as we can see from the computing time for this
% example:
total_time_in_seconds = toc

%%
% For faster work, the first author has been
% exploring electron problems in Julia using Julia's ode45
% command, which runs ten times faster than Matlab's ode45.

