%% Bivariate rootfinding for pattern formation
% Nick Trefethen, October 2019

%%
% (Chebfun example roots/Subramanian.m)
% [Tags: #ROOTS]

%%
% Priya Subramanian at Oxford is interested in the 
% patterns that arise in nonlinear reaction-diffusion PDEs such
% as the Swift-Hohenberg equation and its relatives [2,3,4].  A particular
% interest of hers is cases where the patterns may have quasicrystalline
% structure.  

%%
% In her analysis, systems of polynomial equations arise
% whose roots need to be computed.  She uses the Bertini software
% for this, based on a homotopy method [1].  Sometimes there are just
% two independent variables, however, and this gives nice problems for
% Chebfun2 (chapter 14 of the _Chebfun Guide_).

%%
% Here is one of her examples.  We have parameters $Q$, $\mu$, and $\nu$, for
% which reasonable parameters are these:
Q = 1;
mu = 0.1;
nu = 0.1;

%%
% The independent variables are called $z$ and $w$, and 
% here are the two cubic polynomials of interest:
tic
dom = [-.3 .3 -.15 .15];
z = chebfun2(@(z,w) z, dom);
w = chebfun2(@(z,w) w, dom);
p = mu*z + 2*Q*w.^2 + 4*Q*w.*z - 6*w.^3 - 42*w.^2.*z - 18*w.*z.^2 - 27*z.^3;
q = nu*w + 4*Q*w.*z + 2*Q*z.^2 - 27*w.^3 - 18*w.^2.*z - 42*w.*z.^2 - 6*z.^3;

%%
% Let's plot the zero curves in blue for $p$ and red for $q$,
% with black dots for the common roots:
MS = 'markersize'; LW = 'linewidth';
plot(roots(p),'b',LW,2), hold on, grid on
plot(roots(q),'r',LW,2)
r = roots(p,q)
plot(r(:,1),r(:,2),'.k',MS,20), axis equal, hold off
xlabel z, ylabel w

%%
% For comparison, here is what we get if we negate $\mu$:
mu = -0.1;
p = mu*z + 2*Q*w.^2 + 4*Q*w.*z - 6*w.^3 - 42*w.^2.*z - 18*w.*z.^2 - 27*z.^3;
plot(roots(p),'b',LW,2), hold on, grid on
plot(roots(q),'r',LW,2)
r = roots(p,q)
MS = 'markersize';
plot(r(:,1),r(:,2),'.k',MS,20), axis equal, hold off
xlabel z, ylabel w

%%
Time_for_this_example = toc

%% 
%
% [1] D. J. Bates, J. D. Hauenstein, A. J. Sommese, and
% C. W. Wampler, _Numerically Solving Polynomial Systems with
% Bertini_, SIAM, 2013.
%
% [2] H. Montanelli, Swift-Hohenberg equation in 2D,
% |www.chebfun.org/examples/pde/SwiftHohenberg.html|.
%
% [3] P. Subramanian, A. J. Archer, E. Knobloch, and A. M. Rucklidge,
% Three-dimensional phase field quasicrystals, 
% _Physical Review Letters_ 117:1075501, 2016.
%
% [4] P. Subramanian and A. M. Rucklidge, Mode interactions and
% complex spatial patterns II. Quasicrystals, in preparation.
