%% CF approximation 30 years ago
% Nick Trefethen, July 2014

%%
% (Chebfun example approx/CF30.m)
% [Tags: #Caratheodory-Fejer, #CF]

%%
% In 1986 I published a 4-page proceedings paper called
% "Matlab programs for CF approximation" [1].  It is possible
% that this was the first paper ever published that contained
% Matlab programs.  It presents two programs: CF, for complex
% polynomial or rational CF approximation on the unit disk,
% and RCF, for real polynomial or rational CF approximation on
% the unit interval.

%%
% The real example of CF approximation given in [1]
% is type $(1,1)$ approximation to $(1.2-x)^{1/2}$, for
% which the result is printed as
% 
% $$ r(x) = {1.10417 - 0.77197x \over 1 - 0.27354 x} . $$
%
% In Chebfun, we can perform the calculation like this:
x = chebfun('x');
f = sqrt(1.2-x);
[p,q] = cf(f,1,1);

%%
% The coefficients match nicely:
chebcoeffs(p)
chebcoeffs(q)

%%
% How long does it take?  In [1], the computation took
% "about 40 seconds on an IBM PC/AT".  In Chebfun today, here's
% the timing:
tic, [p,q] = cf(f,1,1); toc

%%
% As far as the eye can tell, the error curve equioscillates between
% 4 extremes:
LW = 'linewidth'; FS = 'fontsize';
errfun = f-p./q; error = norm(errfun,inf);
plot(errfun,LW,1.6)
grid on, ylim(0.02*[-1 1]), hold on
plot([-1 1],error*[+1 +1],'--k',LW,1.6)
plot([-1 1],error*[-1 -1],'--k',LW,1.6), hold off
title(['type (1,1) CF approximation:  error = ' num2str(error)],FS,12)

%%
% Let's compare with the best approximation.  It takes 
% a bit longer to compute:
tic, [p,q] = remez(f,1,1); toc

%%
% Now in principle we have perfect equioscillation,
% but the error is only very slightly smaller:
errfun = f-p./q; error = norm(errfun,inf);
plot(errfun,'m',LW,1.6)
grid on, ylim(0.02*[-1 1]), hold on
plot([-1 1],error*[+1 +1],'--k',LW,1.6)
plot([-1 1],error*[-1 -1],'--k',LW,1.6), hold off
title(['type (1,1) best approximation:  error = ' num2str(error)],FS,12)

%% References
%
% 1. L. N. Trefethen, "Matlab programs for CF approximation,"
%    in C. K. Chui, L. L. Schumaker, and J. D. Ward, eds.,
%    _Approximation Theory V_, Academic Press, 1986, pp. 599-602.
