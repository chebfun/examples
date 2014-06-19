%% Matched asymptotics and boundary layers
% Nick Trefethen, December 2010

%%
% (Chebfun example ode-linear/MatchedAsymp.m)
% [Tags: #linearODE, #asymptotics]

%%
% A powerful technique for problems with large or small parameters is the
% method of matched asymptotics, where approximate solutions accurate in one
% region of the problem domain are matched to different approximate solutions
% accurate in another region.  This subject is discussed in many books by
% authors including Bender & Orszag, Fowler, Howison, Lagerstrom, Nayfeh, and
% van Dyke.

%%
% For example, consider the linear boundary-value problem
%
% $$ -\varepsilon y'' + (2-x^2)y = 1, \qquad    y(-1) = y(1) = 0 $$
%
% with $\varepsilon \ll 1$. In Chebfun, we can set up the problem conveniently
% with a couple of anonymous functions:
d = [-1,1];
x = chebfun('x',d);
L = @(eps) chebop(@(x,u)-eps*diff(u,2) + (2-x.^2).*u,d,0,0);
y = @(eps) L(eps)\1;

%%
% Here are the solutions for three values of $\varepsilon$:
LW = 'linewidth'; FS = 'fontsize';
figure, tic
for j = 1:4
   ep = 10^(-j);
   yep = y(ep);
   subplot(2,2,j), plot(yep,LW,1.6), hold on
   grid on, axis([-1.05 1.05 0 1])
   title(sprintf('eps = %4.1e     npts = %d',ep,length(yep)),FS,8)
end
toc

%%
% It is clear almost at a glance what form the solution is taking as
% $\varepsilon \to 0$. Away from $\pm 1$, the $\varepsilon y''$ term is
% negligible and the solution is approximately that of the rest of the
% equation,
%
% $$ y_{\mbox{interior}} = {1\over 2-x^2}. $$
%
% Near $\pm 1$, on the other hand, $\varepsilon y''$ becomes significant as
% the solution quickly bends down to meet the boundary condition.

%%
% In matched asymptotics the solution away from the boundary layers is called
% the _outer solution_.  Here we have two boundary layers, each of which has
% an _inner solution_. To analyze the right boundary layer, for example, we
% make the approximation $x=1$.  This gives a constant coefficient second-
% order equation, with an exponentially growing solution
% $\exp(x\varepsilon^{-1/2} )$ and an exponentially decaying solution
% $\exp(-x\varepsilon^{-1/2})$. One of our two free parameters is used up by
% the fact that only the first of these is appropriate at the right boundary.
% The other parameter is used to satisfy the boundary condition, giving
%
% $$ y_{\mbox{right}} = 1 - \exp(\varepsilon^{-1/2}(x-1)). $$
%
% Similarly at the left boundary we have
%
% $$ y_{\mbox{left}} = 1 - \exp(\varepsilon^{-1/2}(-x-1)). $$

%%
% The three solutions can be combined to give an asymptotic model valid
% throughout $[-1,1]$:
model = @(eps) 1./(2-x.^2) - exp((x-1)/sqrt(eps)) - exp((-x-1)/sqrt(eps));

%%
% Let us superimpose the prediction of this model, a dashed red line, on the
% plots as before
for j = 1:4
   ep = 10^(-j);
   subplot(2,2,j)
   meps = model(ep);
   plot(meps,'--r',LW,1.6)
   grid on, axis([-1.05 1.05 0 1])
end

%%
% It is interesting to plot and measure the differences between the true
% solution and the model:
for j = 1:4
   ep = 10^(-j);
   subplot(2,2,j)
   yep = y(ep);
   meps = model(ep);
   hold off, plot(meps-yep,'m',LW,1.6)
   grid on, xlim([-1.05 1.05])
   err = norm(yep-meps,inf);
   title(sprintf('eps = %4.1e     err = %5.2e',ep,err),FS,8)
end

%%
% These plots reveal global convergence at a rate $O(\varepsilon^{1/2})$ as
% $\varepsilon\to 0$, with the maximal error being attained in a matching
% region near the boundaries of width $O(\varepsilon^{1/2})$. In the interior
% the accuracy is higher, $O(\varepsilon)$.

%%
% Matched asymptotics is a highly developed field and has been applied to
% linear and nonlinear problems of all kinds. A linear problem with a variable
% coefficient may have interior as well as boundary layers, and for a
% nonlinear problem there may be interior layers at arbitrary locations.
