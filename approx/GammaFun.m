%% The gamma function and its poles
% Nick Hale, December 2009

%%
% (Chebfun example approx/GammaFun.m)
% [Tags: #gamma, #exponents, #pole, #squareroot]

%%
% This script displays some of the features introduced in version 3 for
% unbounded functions by exploring the gamma function $\Gamma(x)$ on the
% interval $[-4,4]$.

%%
% The gamma function has simple poles at the negative integers
% and zero. Chebfun can determine the locations and orders of these poles
% if it is called with the `'blowup'` and `'splitting'` flags on. The
% `'exponents'` field of the output indicates that each pole is simple, that
% is, it has a singularity of type $x^{-1}$.
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
gam = chebfun('gamma(x)',[-4 4],'blowup','on','splitting','on')
plot(gam,'m',LW,1.6), hold on
title('Gamma function',FS,12)

%%
% Alternatively, and always a better idea when the information
% is available, one can instruct Chebfun what poles to put where:
gam = chebfun('gamma(x)',[-4 -3 -2 -1 0 4],'exps',[-1 -1 -1 -1 -1 0])
plot(gam,LW,1.6), hold on
title('Gamma function again',FS,12)

%%
% We can now treat $\Gamma(x)$ like any other chebfun. 
% For example, we can: 
%
% (1) Find its reciprocal $1/\Gamma(x)$:
gam_i = 1./gam;

%%
% (2) Compute the square root $|\Gamma(x)|^{1/2}$:
absgam = abs(gam);
sqrtgam = real(sqrt(absgam));

%%
% (3) Plot these functions:
plot(gam_i,'r', sqrtgam,'-g',LW,1.6)
legend('\Gamma(x)', '1/\Gamma(x)', 'sqrt(|\Gamma(x)|)',...
   'location','southeast')
title('Various related functions',FS,12)

%%
% (4) Plot the critical points:
[y r] = minandmax(gam,'local');
[yi ri] = minandmax(gam_i,'local');
[ys rs] = minandmax(sqrtgam,'local');

plot(r,gam(r),'.k',ri,gam_i(ri),'.k', ...
    rs,sqrtgam(rs),'.k',MS,18,LW,1.6), hold off
title('Gamma function on [-4,4] and its critical points',FS,12)

%%
% (5) Compute some integrals:
sum(gam)
sum(absgam)
sum(sqrtgam)

%%
% Do you understand why these results come out not-a-number, 
% infinite, and finite?
