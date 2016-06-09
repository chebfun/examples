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
N.lbc = [1i; 1];
tic, z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off

%%
% For this computation we have used complex arithmetic for
% convenience, with the nucleus at the origin.  Because of the
% symmetry, only one particle needs to be tracked, so we
% have a scalar complex nonlinear second-order ODE initial-value
% problem:
% $$ z'' = {-2z\over |z|^3} + {1\over z-\overline{z}}. $$
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
% a certain region in the plane.

%% 2. Energy
% The kinetic energy of this motion is $|z'|^2$, and the
% potential energy is $-4/|z| + 1/2\hbox{Im} z$.  Thus the total energy is
% $$ E = |z'|^2 -{4\over |z|} + {1\over 2\hbox{Im} z}, $$
% and this quantity is conserved.  For our initial value
% $V=1$, the energy is $E = -2.5$.  (We do not verify
% this by a Chebfun computation, which would be very slow because
% of the near-singularities when $\hbox{Im}z$ is near zero.)

%%
% Note that a particle at $z=\infty$ with zero velocity
% has energy $0$, and our initial condition will have energy $0$
% with this initial velocity:
% $$ V_{\hbox{crit}} = \sqrt{3.5} \approx 1.87. $$
% Sure enough, one finds that for $V>V_{\hbox{crit}}$, the
% trajectory flies off to infinity.  The reader may find it
% interesting to explore "Pluto" trajectories just below this limit
% starting from values such as $V=1.85$ or $1.86$.

%% 3. Periodic trajectories
% For certain special initial velocities, the orbits are periodic.
% The simplest one, corresponding to $V\approx 1.446$, has
% the electrons simply swinging back and forth:
N.domain = [0 20];
N.lbc = [1i; 1.446];
z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off

%%
% A plot of $x$ values confirms the periodicity:
plot(x), xlabel t, ylabel x

%%
% Here is an estimate of the period:
r = roots(x-.9*max(x)); r = r(deriv(x,r)>0); Period = r(2)-r(1)

%%
% Here is another periodic solution, with an estimate of the period:
N.lbc = [1i; 0.783];
z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off
r = roots(x-.9*max(x)); r = r(deriv(x,r)>0); Period = r(2)-r(1)

%%
% And here is another:
N.lbc = [1i; 1.17745];
z = N\0; x = real(z); y = imag(z);
plot(0,0,'.k',MS,8), hold on
plot(x,y,x,-y,LW,lw), axis(1.2*[-1 1 -1 1]), axis square, hold off
r = roots(y-.99999); Period = mean(r(end-1:end))

%%
% Readers trying these computations on their own may enjoy
% experimenting with the command |comet(z)|.

%% 4. Computing time
% Chebfun is wonderfully convenient, but it is not especially
% fast, as we can see from the computing time for this
% example:
total_time_in_seconds = toc

%%
% For faster work, the first author has been
% exploring electron problems in Julia using Julia's ode45
% command.  In our experiments, this runs more than ten times
% faster than Matlab's ode45.

