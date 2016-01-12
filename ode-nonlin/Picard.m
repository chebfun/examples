%% Picard iteration for ODE existence proof
% Nick Trefethen, January 2016

%%
% (Chebfun example ode-nonlin/Picard.m)
% [Tags: #cumsum, #ode, #Picard]

%%
% In the theory of ODEs, there is a basic theorem of existence
% and uniqueness that is the foundation for everything (see, e.g.,
% [1]).

%%
% *Theorem.*
% _If $f$ is continuous with respect to $t$ and Lipschitz
% continuous with respect to $u$, the first-order IVP_
% $$ u'= f(t,u), \quad t\in [0,d], ~~ \quad u(0) = u_0 $$
% _has a unique solution._

%%
% This theorem applies to systems as well as scalars, and
% since a higher-order equation can be written as a system of
% first-order equations, it covers higher-order ODEs too.

%%
% The standard proof of the theorem is based on Picard (=
% Picard-Lindelof) iteration,
% and it can be illustrated in Chebfun using the
% |cumsum| command.  The proof is based on noting that the ODE
% is equivalent to
% $$ u(t) = u_0 + \int_0^t f(s,u(s)) ds . $$
% In the Picard iteration, we successively consider
%
% $$ u^{(0)} = u_0, ~~~
% u^{(1)} = u_0 + \int_0^t u^{(0)}(s) ds,  ~~~
% u^{(2)} = u_0 + \int_0^t u^{(1)}(s) ds, $$
%
% and so on.  One can prove with the given assumptions that
% this process converges to a unique solution of (1).

%%
% Let us see the iteration in action for the problem
% $$ u' = \sin(u) + \sin(t), \quad t\in [0,8], ~~ u(0) = 1, $$ 
% whose solution we can compute like this:
d = 8; t = chebfun('t',[0 d]); u0 = 1;
L = chebop(0,d); L.op = @(t,u) diff(u) - sin(u); L.lbc = u0;
uexact = L\sin(t);

%%
% This first plot shows iterates $k = 0,\dots,4$, with
% the exact solution in red.
u = u0 + 0*t;
f = @(u,t) sin(u) + sin(t);
LW = 'linewidth'; FS = 'fontsize'; IN = 'interpret'; LT = 'latex';
hold off
ss = @(k) ['$k = ' int2str(k) '$'];
for k = 0:4
  plot(u,'b',LW,1), hold on, ylim([-3 10])
  text(1.015*d,u(end),ss(k),IN,LT)
  u = u0 + cumsum(f(u,t));
end
plot(uexact,'r',LW,1), xlabel('t',FS,10), ylabel('u',FS,10)
title('Picard iterates $k = 0,\dots,4$',FS,12,IN,LT)

%%
% A second plot shows $k = 5,\dots,9$.
hold off
for k = 5:9
  plot(u,'b',LW,1), hold on, ylim([0 7])
  text(1.015*d,u(end),ss(k),IN,LT)
  u = u0 + cumsum(f(u,t));
end
plot(uexact,'r',LW,1), xlabel('t',FS,10), ylabel('u',FS,10)
title('Picard iterates $k = 5,\dots,9$',FS,12,IN,LT)

%%
% A third plot shows $k = 10,\dots ,14$.
hold off
for k = 10:14
  plot(u,'b',LW,1), hold on, ylim([1 6])
  text(1.015*d,u(end),ss(k),IN,LT)
  u = u0 + cumsum(f(u,t));
end
plot(uexact,'r',LW,1), xlabel('t',FS,10), ylabel('u',FS,10)
title('Picard iterates $k = 10,\dots,14$',FS,12,IN,LT)

%%
% These plots show vividly the kind of convergence one
% can expect from a Picard iteration: starting at the
% initial condition, sweeping slowly across the domain.
% There is a numerical method based on this idea, called
% _waveform relaxation_, but one can see immediately from
% the pictures that it is unlikely to be efficient when
% carried out over long time intervals.  Instead, standard
% numerical methods just march once rather than many times
% from left to right,
% but they march with a small discrete time step and a 
% discrete formula of higher order.

%%
% To see the convergence in a quantitative fashion, it is interesting
% plot the errors of iterates $0,\dots, 4$ as a function of $t$ 
% on a loglog plot.  The zeroth iterate has accuracy $O(t)$, the
% first has accuracy $O(t^2)$, and so on:

u = u0 + 0*t;
ss = @(k) ['$k = ' int2str(k) '$'];
tt = logspace(-2,log10(8),600); hold off
for k = 0:4
  errtt = abs(u(tt)-uexact(tt));
  loglog(tt,errtt,'k',LW,.7), hold on
  text(8.7,errtt(1),ss(k),IN,LT)
  u = u0 + cumsum(f(u,t));
end
xlabel('t',FS,10), ylabel('error',FS,10)
axis([1e-2 8 1e-16 1e3])
title('Errors of iterates $0,\dots,4$',FS,12,IN,LT)

%%
% Reference:
%
% [1] E. A. Coddington and N. Levinson, 
% _Theory of Ordinary Differential Equations_, McGraw-Hill, 1955.
