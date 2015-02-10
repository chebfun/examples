%% Location and stability of periodic orbits in the Van der Pol equation
% Toby Driscoll and Hrothgar, April 2014

%%
% (Chebfun Example ode-nonlin/VdPOrbit)
% [Tags: #ODE, #periodic, #stability]

%%
plotargs = {'linewidth',2,'markersize',8};

%% The Van der Pol equation
% The Van der Pol equation is a classic nonlinear oscillator, originally used to
% model a diode:
%
% \[ u'' - \mu(1-u^2) u' + u = 0, \quad u(0), u'(0) \text{ given}. \]
% 
% We can see that the equation has an attractive periodic orbit for any positive
% value of $\mu$. 

mu = 1;
f = @(t,u) diff(u,2) - mu*(1-u.^2).*diff(u) + u;
N = chebop(f,[0,18]);
N.bc = @(t,u) [ u(0)-1, feval(diff(u),0) ];

%%
u = N\0;
plot(u,diff(u),plotargs{:})
set(gcf,'defaulttextinterpreter','tex')
xlabel('u'),  ylabel('u''')
title('Attraction to the limit cycle')

%%
% To find the limit cycle, we could simply solve for a long time. But
% there is a more direct route. The trick is that the period of the cycle is not
% known in advance. In order to pose a problem on a fixed domain, we introduce a
% new time variable $\tau = t/T$, where $T$ is the unknown period, so that
% $\tau\in[0,1]$. We add the equation $d T / d \tau = 0$ to close the system. 

ODE = @(tau,u,T) [ diff(u,2)./T.^2 - mu*(1-u.^2).*diff(u)./T + u , diff(T)];
N = chebop( ODE, [0 1] );

%%
% We want the solution to be periodic, which requires continuity in both value
% and the first derivative of $u$. Finally, we also need an "anchor" or phase
% condition, in order to avoid an infinite family of solutions parameterized by
% the starting point at $\tau=0$. 
N.bc = @(tau,u,T) [ u(0)-u(1), feval(diff(u),0)-feval(diff(u),1), u(0)-1 ];

%%
% The system needs to be seeded with a decent initial guess to converge. 
tau = chebfun(@(t) t,[0,1]);
N.init = [ 1+sin(2*pi*tau), chebfun(5,[0,1]) ];

%%
% We'll lighten the accuracy requirements to save computaion time.
cheboppref('restol',1e-5,'deltol',1e-5)
chebfunpref('eps',1e-8)
cheboppref('display','iter')

z = N\0;

%%
u = z(:,1);
period = z(0,2)
hold on, plot(u,diff(u)/period,'r',plotargs{:})

%%
% For larger values of $\mu$, it's helpful to continue from each found solution
% as the starting guess. 
z = {z};
for mu = [2,5,10]
    N.op = @(tau,u,T) [ diff(u,2)./T.^2 - mu*(1-u.^2).*diff(u)./T + u , diff(T)];
    N.init = z{end};
    z{end+1} = N\0;
    u = z{end}(:,1);
    period = z{end}(0,2);
    hold on,  plot(u,diff(u)/period,'k',plotargs{:})
    fprintf(' mu = %i, period = %.6f\n\n',mu,period)
end

%%
% For this value of $\mu$, the solution has some pretty steep gradients, which
% makes a global solution less than ideal.

%% Stability and the monodromy matrix
% Let's return to a gentler periodic solution for a smaller $\mu$. We are going
% to regard the period $T$ as fixed now, and represent our solutions on the
% natural domain. 
mu = 2;  z = z{2};
T = z(0,2);
u = chebfun( @(t) z(t/T,1), [0 T] );

%%
% From the dynamics, we would think that this solution is asymptotically stable,
% meaning that perturbations to it will return to the orbit. But there is a
% different tool for testing that idea more thoroughly: the _monodromy matrix_ [1].

%%
% We begin by rewriting the solution in the form of a first-order system.
% Define $v=u'$, so that the derivative of $[u;v]$ is 
F = @(u,v) [ v, mu*(1-u.^2).*v - u ];

%%
% Now we will linearize the dynamics around the periodic solution using
% Chebfun's automatic differentiation. This requires a subtle trick. From one
% point of view, $v$ is derived from $u$, so derivatives of $F$ all are as well.
% But we want to tell Chebfun instead to regard $v$ as an independent variable
% in the first order system. We use JACRESET to do this.
v = diff(u);
v = jacreset(v);

%%
% Now we can easily compute the Jacobian (or Frechet derivative) of the system
% around the solution. We will find it convenient to store the columns
% separately. 
r = F(u,v);
Ju = diff(r,u);
Jv = diff(r,v);

%%
% The monodromy matrix is defined as the matrix solution of $Y'(t)=J(t)Y(t)$,
% $Y(0)=I$, evaluated at the final time $t=T$. 
linode =  @(t,y) diff(y) - (Ju*y(:,1) + Jv*y(:,2));
A = chebop(linode, [0,T] );
A.init = [ 0*u, 0*v ];

%%
% Use the first column of the identity as the initial condition:
A.bc = @(t,y) [y(0,1)-1, y(0,2)];
Y1 = A\0;

%%
% Now, use the second column of the identity as the initial condition:
A.bc = @(t,y) [y(0,1), y(0,2)-1];
Y2 = A\0;

subplot(2,1,1)
plot(Y1,plotargs{:}),  title('Solutions in first column')
subplot(2,1,2)
plot(Y2,plotargs{:}),  title('Solutions in second column')

%% 
% We assemble the monodromy matrix and look at its eigenvalues.
M = [ Y1(T,:)', Y2(T,:)' ]

%%
format long
eig(M)

%%
% Theoretically, one eigenvalue of $M$ should be exactly one (recall that we
% only asked for about 5 digits of accuracy in the periodic solution), and the
% others tell the story of stability. The second eigenvalue being so small
% indicates that the limit cycle is very stable, even though the linearized
% departures from the cycle can grow quite a bit during one loop. 

%% Negative $\mu$: an unstable cycle
% If $\mu < 0$ in the Van der Pol equation, we can still find a periodic
% solution. 

mu = -1;
cheboppref factory
chebfunpref factory
ODE = @(tau,u,T) [ diff(u,2)./T.^2 - mu*(1-u.^2).*diff(u)./T + u , diff(T)];
N = chebop( ODE, [0 1] );
N.bc = @(tau,u,T) [ u(0)-u(1), feval(diff(u),0)-feval(diff(u),1), u(0)-1 ];
N.init = [ 1+sin(2*pi*tau), chebfun(5,[0,1]) ];

%%
z = N\0;
u = z(:,1);
T = z(0,2)
clf, plot(u,diff(u)/period,plotargs{:})
xlabel('u'),  ylabel('u''')
title('Periodic orbit for negative \mu')


%%
% We again linearize to prepare for computing the monodromy matrix.
u = chebfun( @(t) z(t/T,1), [0 T] );
v = jacreset( diff(u) );
r = [ v, mu*(1-u.^2).*v - u ];
Ju = diff(r,u);  Jv = diff(r,v);
linode =  @(t,y) diff(y) - (Ju*y(:,1) + Jv*y(:,2));
A = chebop(linode, [0,T] );
A.init = [ 0*u, 0*v ];

%%
% Now, one of the eigenvalues is much larger than unity:
A.bc = @(t,y) [y(0,1)-1, y(0,2)];
Y1 = A\0;

A.bc = @(t,y) [y(0,1), y(0,2)-1];
Y2 = A\0;

M = [ Y1(T,:)', Y2(T,:)' ]

%%
eig(M)

%%
% We conclude that the periodic solution in this case is asymptotically
% unstable. In this simple case, we could find it by solving Van der Pol
% backwards in time, which is equivalent to reversing the sign of $\mu$. 


%% Reference
% 1. R. Seydel, "New methods for calculating the stability of periodic
%    solutions." Comp. Math. Appl. 14(7), 505-510, 1987.
