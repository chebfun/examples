%% Pitchfork bifurcation triggered by noise
% Nick Trefethen, May 2017

%%
% (Chebfun example ode-random/Bifurcation.m)

%%
% The second-order ODE
% $$ y'' = 2\kern .3pt cy - 4y^3 - 0.2\kern .3pt y' $$
% has steady states $y=0$, $y = \sqrt{c/2}$ and $y = \sqrt{-c/2}$.
% For $c<0$, only the first of these is real, and it is a stable
% steady state.  As $c$ passes through zero to values $c>0$, 
% the other two steady states emerge, and now $y=0$ is unsteady
% and the other two are steady.  This is a pitchfork bifurcation.

%%
% One way to see the bifurcation in a single ODE computation
% is to consider a time-dependent problem with a coefficient $c(t)$
% that slowly increases through zero, like this:
% $$ y'' = 2\kern .3pt c(t)y - 4y^3 - 0.2\kern .3pt y' ,
% ~~~ c(t) = -1+t/300, ~~ t \in [0, 600].  $$
% If the initial condition is $y(0) = 0$, then the solution is
% $y(t)=0$ for all $t$ --- the solution never notices the instability.

%%
% However, suppose we add noise in the form of a random function
% of small amplitude, like this:
% $$ y'' = 2\kern .3pt c(t)y - 4y^3 - 0.2\kern .3pt y' + 0.003f(t),
% ~~~ c(t) = -1+t/300, ~~ t \in [0, 600].  $$
% Now the solution will, at random, deviate to one or the
% other branch of the pitchfork.  Here we show three solutions:
% two with noise, and one without (the dashed middle line).
% (We tried |rng(0) ... rng(6)| before finally finding that
% |rng(7)| gave a picture with both forks!)
tic
rng(7), lambda = 2;
N = chebop(0,600); N.lbc = [0;0];
N.op = @(t,y) diff(y,2) - 2*(-1+t/300)*y + 4*y^3 + .2*diff(y);
FS = 'fontsize'; LW = 'linewidth';
y1 = N\0; plot(y1,'--k',LW,2), hold on
f1 = 0.003*randnfun([0 600],lambda,'norm');
y2 = N\f1; plot(y2,'b',LW,2),
f2 = 0.003*randnfun([0 600],lambda,'norm');
y3 = N\f2; plot(y3,'r',LW,2), hold off
xlabel('t',FS,32), ylabel('y',FS,32)
title('Pitchfork with damping',FS,32)
axis([0 600 -.8 .8]), grid on

%%
% Note that the ODE has a damping term involving $y'$.
% If we remove the damping, we get a similar picture but with
% large oscillations:
N.op = @(t,y) diff(y,2) - 2*(-1+t/300)*y + 4*y^3;
plot(y1,'--k',LW,2), hold on
y2 = N\(-f1); plot(y2,'b',LW,2),
y3 = N\f2; plot(y3,'r',LW,2), hold off
xlabel('t',FS,32), ylabel('y',FS,32)
title('Pitchfork without damping',FS,32)
axis([0 600 -.8 .8]), grid on

%%
% When we first tried this, both branches went positive.
% So we flipped the sign on one of them.

%%
% This example has some close links with the
% example "Phase-locking in a Duffing-type equation", though
% that involves a first-order ODE.

%%
total_time_in_seconds = toc
